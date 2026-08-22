# Legacy Python prototypes

These files are the original standalone Python models the C++ library (and its
[Python API](../../README.md#python-api)) superseded. They are **kept for historical reference
only** and are not maintained, not tested, and not guaranteed to run (several have bugs: missing
`return`s, undefined attributes, mismatched CoolProp call signatures).

| Legacy file | C++ / Python API equivalent |
|---|---|
| `reactor.py`, `reactor_model_B.py` | `astara::reactor::Reactor` / `astara.Reactor` |
| `Pressurizer.py` | `astara::pressurizer::Pressurizer` / `astara.Pressurizer` |
| `PrimaryCoolantPump.py`, `reactor_primary_coolant_pump.py` | `astara::pump::HomologousPump` / `astara.HomologousPump` |
| `SteamGenerator.py`, `steamgeneratormodelB.py` | `astara::sg::AliSteamGenerator` / `astara::sg::HelicalCoilSteamGenerator` (`astara.AliSteamGenerator` / `astara.HelicalCoilSteamGenerator`) |
| `CondenserModel.py` | `astara::condenser::Condenser` / `astara.Condenser` |
| `ControlRodSystem.py` | `astara::control::ControlRodSystem` / `astara.ControlRodSystem` |
| `Base.py` | superseded by `astara::core::Integrator` (`rk4Step`) |

`TurbineModel.py` and `turbine_new.py` remain **unported** -- there is no C++ or Python API
equivalent yet. Both files are themselves incomplete/inconsistent prototypes (undefined
coefficients such as `Kch`, several "at random" / "need to calculate" placeholder constants,
inconsistent cross-object wiring between the two files' approaches), so porting them is left as
future work rather than a straight translation.
