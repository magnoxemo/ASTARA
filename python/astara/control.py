"""Control-rod ergonomics on top of the bound `ControlRodSystem`.

`astara.ControlRodSystem.reactivity_rate(mode)` (see
`include/astara/control/ControlRodSystem.hpp`) only computes a rate -- the
caller is responsible for integrating `rate * dt` into the target reactor's
`reactor.reactivity().rho_external` every step, and for remembering which
mode (Insert/Withdraw/Hold) is currently commanded. `Rod` does both, so a
notebook can just call `rod.insert()` / `rod.withdraw()` / `rod.hold()` and
step the plant -- no manual reactivity bookkeeping.

`ControlRodSystem` itself has no notion of rod position -- its rate is
normalised (a fraction of full travel per second), so `Rod` integrates it
into a `position` in [0.0, 1.0] (0.0 = fully inserted, 1.0 = fully
withdrawn, starting at 0.5) and clamps there: once a rod bottoms/tops out,
further motion in that direction stops contributing reactivity, and a
`RuntimeWarning` fires (once per limit reached, not once per step).
"""

from __future__ import annotations

import warnings

from ._astara import ControlRodSystem, ControlRodSystemConfig, RodMode


class Rod(object):
    """A control-rod bank bound to one reactor.

    ```python
    rod = astara.Rod(reactor)
    rod.insert()          # commands insertion; integrated on every time_step()
    ...
    rod.hold()
    rod.position           # 0.0 (fully inserted) .. 1.0 (fully withdrawn)
    ```

    Add a `Rod` to a `Plant` like any other component -- `Plant.step()` calls
    `time_step(dt)` on everything it holds, and `Rod.time_step` integrates the
    commanded rate into the bound reactor's external reactivity, clamped to
    the rod's travel range.
    """

    def __init__(self, reactor, config: ControlRodSystemConfig | None = None):
        self._reactor = reactor
        self._system = ControlRodSystem(config or ControlRodSystemConfig())
        self._mode = RodMode.Hold
        self._position = 0.5
        self._limit_reported: str | None = None

    @property
    def config(self) -> ControlRodSystemConfig:
        return self._system.config()

    @property
    def mode(self) -> RodMode:
        return self._mode

    @property
    def position(self) -> float:
        """Normalised rod-bank position: 0.0 = fully inserted, 1.0 = fully withdrawn."""
        return self._position

    def insert(self) -> None:
        self._mode = RodMode.Insert

    def withdraw(self) -> None:
        self._mode = RodMode.Withdraw

    def hold(self) -> None:
        self._mode = RodMode.Hold

    def time_step(self, dt: float) -> None:
        """Integrate the commanded reactivity rate into the bound reactor.

        Named `time_step` (not e.g. `advance`) so a `Plant`/`Simulation` can
        step a `Rod` with the same call it uses for every other component.
        Motion (and the reactivity it contributes) stops at `position` 0.0
        or 1.0 -- see the module docstring.
        """
        if self._mode is RodMode.Hold:
            return

        cfg = self._system.config()
        position_rate = -cfg.insertion_speed if self._mode is RodMode.Insert else cfg.withdrawal_speed
        target_position = self._position + position_rate * dt
        new_position = min(1.0, max(0.0, target_position))
        moved = new_position - self._position

        if moved != 0.0:
            reactivity = self._reactor.reactivity()
            reactivity.rho_external = reactivity.rho_external + cfg.step_reactivity * moved
            self._position = new_position

        if new_position == target_position:
            self._limit_reported = None
            return

        limit = "fully withdrawn" if new_position >= 1.0 else "fully inserted"
        if self._limit_reported != limit:
            warnings.warn(
                f"Rod {limit} (position={new_position:.3f}); "
                f"{self._mode.name.lower()} has no further effect until reversed",
                RuntimeWarning,
                stacklevel=2,
            )
            self._limit_reported = limit
