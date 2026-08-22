"""Python bindings for the ASTARA PWR primary-loop simulator.

Re-exports the compiled `_astara` nanobind extension module (all physics)
under a clean `astara.*` namespace, plus a pure-Python ergonomics layer on
top of it:

- `Plant` / `Loop` (`astara.network`) -- wire arbitrary numbers of
  components together (e.g. several reactor+pump loops fanned into one
  shared steam generator), instead of being limited to the single fixed
  1-reactor:1-SG:1-pump:1-pressurizer topology of `PrimaryLoop`.
- `Rod` (`astara.control`) -- insert/withdraw/hold a control-rod bank on a
  reactor without manually integrating a reactivity rate every step.
- `Simulation` (`astara.sim`) -- a run loop with scheduled one-off events
  and DataFrame recording, replacing the hand-rolled `for i in
  range(n_steps)` loops every demo used to write.
- `astara.postprocess` -- quick-plot and comparison helpers for the
  DataFrames `Simulation.run()` returns.

The `network`/`control`/`sim` layer is pure orchestration around the
existing bound classes -- no new C++.
"""

from ._astara import *  # noqa: F401,F403
from ._astara import __doc__ as _doc  # noqa: F401

from .network import Plant, Loop, register_adapter  # noqa: F401
from .control import Rod  # noqa: F401
from .sim import Simulation  # noqa: F401
from . import postprocess  # noqa: F401

__all__ = [name for name in dir() if not name.startswith("_")]
