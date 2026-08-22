"""Flexible plant topologies: wire components together without touching C++.

The compiled `astara.PrimaryLoop` couples exactly one reactor, one steam
generator, one pump, and one pressurizer. Real plants aren't always that
shape -- e.g. several small-modular-reactor loops (reactor + pump each)
feeding a single shared steam generator. `Plant` builds that kind of
topology out of the same per-component bindings `PrimaryLoop` itself uses,
generalized to fan-in (many sources -> one destination) and fan-out (one
source -> many destinations).

This is pure-Python orchestration: every step it reads each source
component's outlet quantities, merges them, and writes the result into each
destination's inlet. It carries the same fidelity as `PrimaryLoop::timeStep`
(temperature and mass flow only, no enthalpy/pressure coupling) -- see
`include/astara/primary/PrimaryLoop.hpp` for the equivalent single-loop
version this generalizes.

Mixing rule at a fan-in point (multiple sources -> one destination):
total mass flow = sum of the branch flows; merged temperature = the
mass-flow-weighted average of the branch temperatures (a branch with no
flow reading of its own is weighted equally, i.e. weight 1.0). This is the
standard lumped-node mixing assumption, not a full enthalpy balance.
"""

from __future__ import annotations

from typing import Callable, Optional


class Loop:
    """One reactor + pump (+ optional pressurizer): the repeatable SMR-style unit.

    Bundles the objects that together produce one thermal-hydraulic branch
    (hot-leg temperature from the reactor, mass flow from the pump) and
    accepts the returning cold-leg temperature on the reactor. Use several
    `Loop`s as the sources of a `Plant.connect()` fan-in to a shared steam
    generator.
    """

    def __init__(self, reactor, pump, pressurizer=None, name: Optional[str] = None):
        self.reactor = reactor
        self.pump = pump
        self.pressurizer = pressurizer
        self.name = name

    def time_step(self, dt: float) -> None:
        self.reactor.time_step(dt)
        self.pump.time_step(dt)
        if self.pressurizer is not None:
            self.pressurizer.time_step(dt)

    def get_outlet(self) -> dict:
        return {
            "temperature_K": self.reactor.hot_leg_outlet_temperature_K(),
            "mass_flow_kg_s": self.pump.mass_flow_kg_s(),
        }

    def set_inlet(self, **fields) -> None:
        if "temperature_K" in fields:
            self.reactor.set_cold_leg_inlet_temperature(fields["temperature_K"])


def _reactor_outlet(reactor) -> dict:
    return {"temperature_K": reactor.hot_leg_outlet_temperature_K()}


def _reactor_inlet(reactor, **fields) -> None:
    if "temperature_K" in fields:
        reactor.set_cold_leg_inlet_temperature(fields["temperature_K"])


def _pump_outlet(pump) -> dict:
    return {"mass_flow_kg_s": pump.mass_flow_kg_s()}


def _sg_outlet(sg) -> dict:
    return {"temperature_K": sg.primary_outlet_temperature_K()}


def _sg_inlet(sg, **fields) -> None:
    inputs = sg.inputs()
    if "temperature_K" in fields:
        inputs.primary_inlet_temperature_K = fields["temperature_K"]
    if "mass_flow_kg_s" in fields:
        inputs.primary_mass_flow_kg_s = fields["mass_flow_kg_s"]


# Maps a bound component's *type name* to (get_outlet, set_inlet) functions.
# Keyed by name (rather than the class object) so this doesn't need to import
# the compiled `_astara` extension types directly; entries are looked up by
# each object's `type(obj).__name__` at connect()/step() time.
_ADAPTERS: dict[str, tuple[Callable, Callable]] = {
    "Reactor": (_reactor_outlet, _reactor_inlet),
    "HomologousPump": (_pump_outlet, lambda pump, **_fields: None),
    "AliSteamGenerator": (_sg_outlet, _sg_inlet),
    "HelicalCoilSteamGenerator": (_sg_outlet, _sg_inlet),
    "Loop": (Loop.get_outlet, Loop.set_inlet),
}


def register_adapter(type_name: str, get_outlet: Callable, set_inlet: Callable) -> None:
    """Teach `Plant.connect()`/`step()` how to read/write a new component type."""
    _ADAPTERS[type_name] = (get_outlet, set_inlet)


def _adapter_for(obj):
    type_name = type(obj).__name__
    try:
        return _ADAPTERS[type_name]
    except KeyError as exc:
        raise TypeError(
            f"No network adapter registered for {type_name!r}. "
            "Use astara.network.register_adapter() to add one, or connect a "
            "Loop/Reactor/Pump/SteamGenerator instead."
        ) from exc


def _merge_outlets(sources: list) -> dict:
    outlets = [_adapter_for(src)[0](src) for src in sources]
    merged: dict = {}

    flows = [o["mass_flow_kg_s"] for o in outlets if "mass_flow_kg_s" in o]
    if flows:
        merged["mass_flow_kg_s"] = sum(flows)

    temps = [(o["temperature_K"], o.get("mass_flow_kg_s", 1.0)) for o in outlets if "temperature_K" in o]
    if temps:
        total_weight = sum(w for _, w in temps)
        merged["temperature_K"] = sum(t * w for t, w in temps) / total_weight if total_weight else (
            sum(t for t, _ in temps) / len(temps)
        )

    return merged


class Plant:
    """A named collection of components wired together with `connect()`.

    ```python
    plant = astara.Plant()
    loops = [plant.add(f"loop{i}", astara.Loop(make_reactor(), make_pump()))
             for i in range(4)]
    sg = plant.add("sg", make_sg(props))
    plant.connect(loops, sg)   # fan-in:  mixed T + summed flow -> sg.inputs()
    plant.connect(sg, loops)   # fan-out: sg outlet T -> every loop's cold leg

    for _ in range(n_steps):
        plant.step(dt)
    ```
    """

    def __init__(self):
        self._components: dict[str, object] = {}
        self._edges: list[tuple[list, list]] = []

    def add(self, name: str, component):
        """Register `component` under `name` and return it (for chaining)."""
        self._components[name] = component
        return component

    def __getitem__(self, name: str):
        return self._components[name]

    def __iter__(self):
        return iter(self._components.values())

    def connect(self, source, dest) -> None:
        """Wire `source` -> `dest`. Either side may be a single component or a
        list, giving fan-in (list -> one) or fan-out (one -> list) coupling."""
        sources = source if isinstance(source, (list, tuple)) else [source]
        dests = dest if isinstance(dest, (list, tuple)) else [dest]
        for obj in list(sources) + list(dests):
            _adapter_for(obj)  # validates eagerly, fails fast on typos
        self._edges.append((list(sources), list(dests)))

    def step(self, dt: float) -> None:
        """Advance every component by `dt`, then propagate connections."""
        for component in self._components.values():
            component.time_step(dt)
        for sources, dests in self._edges:
            merged = _merge_outlets(sources)
            for dest in dests:
                _adapter_for(dest)[1](dest, **merged)
