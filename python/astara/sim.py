"""Run loop + event scheduling + recording, in one object.

Every demo script before this module hand-rolled the same three things: a
`for i in range(n_steps)` loop, an `if not fired and t >= t_trigger:` guard
for one-off events (pump trips, rod moves), and a CSV-writing block for
logging. `Simulation` replaces all three with `schedule()`, `record()`, and
`run()`, returning a `pandas.DataFrame` (or a plain list of dicts if
`pandas` isn't installed) ready for `astara.postprocess` or a notebook's own
plotting.
"""

from __future__ import annotations

from typing import Callable


class Simulation:
    """Drives a `Plant` (or any object with `.step(dt)`/`.time_step(dt)`).

    ```python
    sim = astara.Simulation(plant, dt=1e-3)
    sim.schedule(t=20.0, action=lambda: pump.set_input_power_W(0.0))
    sim.record(lambda: {"power": reactor.state().kinetics.power}, every=0.1)
    df = sim.run(t_end=60.0)
    ```
    """

    def __init__(self, plant, dt: float):
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        self.plant = plant
        self.dt = dt
        self._scheduled: list[dict] = []
        self._recorders: list[dict] = []

    def schedule(self, t: float, action: Callable[[], None]) -> None:
        """Call `action()` exactly once, the first time simulated time reaches `t`."""
        self._scheduled.append({"t": t, "action": action, "fired": False})

    def record(self, fn: Callable[[], dict] | dict[str, Callable[[], float]], every: float) -> None:
        """Sample `fn` every `every` seconds and stash the result as one row.

        `fn` is either a callable returning a flat dict, or a dict mapping
        each column name straight to a zero-arg getter, e.g.:

        ```python
        sim.record({"T_hot": reactor.hot_leg_outlet_temperature_K, "flow": pump.mass_flow_kg_s}, every=0.1)
        ```

        -- pass whichever is more convenient; both produce identical rows,
        so a per-loop dict comprehension replaces a hand-written wrapper
        function. Multiple `record()` calls are independent and may use
        different cadences. The row does not need a "t_s" key -- it's added
        automatically.
        """
        if every <= 0.0:
            raise ValueError("every must be positive")
        if isinstance(fn, dict):
            fields = fn
            fn = lambda: {name: getter() for name, getter in fields.items()}
        self._recorders.append({"fn": fn, "every": every, "next_t": 0.0, "rows": []})

    def _step_plant(self, dt: float) -> None:
        step = getattr(self.plant, "step", None)
        if step is not None:
            step(dt)
        else:
            self.plant.time_step(dt)

    def run(self, t_end: float):
        """Advance from t=0 to `t_end`, firing scheduled actions/recorders.

        Returns a single `pandas.DataFrame` if exactly one `record()` was
        registered, a list of DataFrames (one per `record()` call, in
        registration order) if there were several, or `None` if none were.
        Falls back to plain lists of dicts wherever `pandas` isn't installed.
        """
        dt = self.dt
        n_steps = int(round(t_end / dt))

        for i in range(n_steps + 1):
            t = i * dt

            for event in self._scheduled:
                if not event["fired"] and t >= event["t"]:
                    event["action"]()
                    event["fired"] = True

            for rec in self._recorders:
                if t + 0.5 * dt >= rec["next_t"]:
                    row = dict(rec["fn"]())
                    row.setdefault("t_s", t)
                    rec["rows"].append(row)
                    rec["next_t"] += rec["every"]

            if i < n_steps:
                self._step_plant(dt)

        frames = [_to_frame(rec["rows"]) for rec in self._recorders]
        if not frames:
            return None
        return frames[0] if len(frames) == 1 else frames


def _to_frame(rows: list[dict]):
    try:
        import pandas as pd
    except ImportError:
        return rows
    return pd.DataFrame(rows)
