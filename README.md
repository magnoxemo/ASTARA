# ASTARA

ASTARA C++ dynamic simulation library for the primary loop of a pressurized-water reactor (PWR), 
with extensions for small modular reactors (SMRs) using helical-coil once-through 
steam generators. Also, the models are generalized enough to support other designs as
well. 

The over all design followed an OOP design and should be easily customizable. 


![CI](https://github.com/du-ards/ASTARA/actions/workflows/workflow.yml/badge.svg?branch=main)
![Docker](https://github.com/du-ards/ASTARA/actions/workflows/docker_build.yml/badge.svg?branch=main)


# Dependencies

For the thermal fluid property lookup ASTARA relies on a third party lib called 
[IF97](https://github.com/coolprop/IF97). The code is tested gcc-13.3.0 
version in Ubuntu 24.04 LTS and Debian version 13.0. For the build system it relies entirely 
on CMake (VERSION 3.15 or later ). For Windows support please 
reach out to any of the maintainers. 

# Setup

If you have a linux machine then you are already ready to go. 
Just copy the commands and run one by one.

```
git clone https://github.com/du-ards/ASTARA.git
cd ASTARA
git submodule update --init --recursive 

mkdir -p build 
cd build 
cmake ..
make -j $nproc
```
# Adding as a CMake project

We recommend that add this repo as a submodule and then add the repositiory as a subdirectory 
in your cmake build system. 

# Python API

ASTARA also ships a [nanobind](https://github.com/wjakob/nanobind)-based Python API that binds
the C++ library directly -- no separate Python physics implementation, all computation stays in
C++. It supersedes the standalone prototypes in `legacy/python/` (see
[`legacy/python/README.md`](legacy/python/README.md) for the mapping).

Install with pip (builds the extension via scikit-build-core + CMake under the hood):

```
git clone --recurse-submodules https://github.com/du-ards/ASTARA.git
cd ASTARA
pip install .              # core bindings only
pip install .[analysis]    # + pandas/matplotlib, needed for Simulation recording/postprocess
```

Or, for iterating on the bindings themselves, build straight from CMake and point `PYTHONPATH`
at the build tree instead of installing:

```
git submodule update --init --recursive   # pulls in externals/nanobind too
mkdir -p build && cd build
cmake .. -DASTARA_BUILD_PYTHON=ON
make -j $nproc

PYTHONPATH=$(pwd)/python python3 -c "import astara"
```

Every C++ module is bound: `core`, `props` (`IF97Water`, `LinearizedWater`), `reactor`, `pump`,
`pressurizer`, `sg` (`AliSteamGenerator`, `HelicalCoilSteamGenerator`), `control` (including
`ControlRodSystem`), and `primary.PrimaryLoop`. Run the Python smoke-test suite with
`ctest -R python_smoke` from the build directory (needs `pytest`).

## Composing plants, rods, and post-processing from Python

On top of that 1:1 binding layer, `astara` adds pure-Python ergonomics aimed at running whole
scenarios from a notebook with the least code possible:

- **`astara.Plant` / `astara.Loop`** -- wire arbitrary numbers of components together instead of
  being limited to `PrimaryLoop`'s fixed one-reactor-one-SG-one-pump-one-pressurizer shape. e.g.
  several SMR reactor+pump loops fanned into one shared steam generator (see
  `python/demo/four_smr_one_sg.py`).
- **`astara.Rod`** -- `rod.insert()` / `rod.withdraw()` / `rod.hold()` on a reactor; the reactivity
  rate is integrated automatically every step, no manual bookkeeping.
- **`astara.Simulation`** -- a run loop with `schedule(t, action)` for one-off events (pump trips,
  rod moves) and `record(fn, every)` for sampling, returning a `pandas.DataFrame` from `run()`.
- **`astara.postprocess`** -- `plot_traces()`, `compare()` (baseline-vs-perturbed), `to_csv()` /
  `to_parquet()` on those DataFrames.

```python
import astara

plant = astara.Plant()
loops = [plant.add(f"loop{i}", astara.Loop(make_reactor(), make_pump())) for i in range(4)]
sg = plant.add("sg", make_sg(props))
plant.connect(loops, sg)   # fan-in:  4 loops -> shared SG inlet (flow-weighted mixing)
plant.connect(sg, loops)   # fan-out: shared SG outlet -> every loop's cold leg

rod = astara.Rod(loops[0].reactor)
sim = astara.Simulation(plant, dt=1.0e-3)
sim.schedule(t=10.0, action=rod.insert)
sim.record(lambda: {"power0": loops[0].reactor.state().kinetics.power}, every=0.1)

df = sim.run(t_end=30.0)
astara.postprocess.plot_traces(df, ["power0"])
```

See `python/demo/four_smr_one_sg.py` for the full runnable version of this example (plus the
`make_reactor`/`make_pump`/`make_sg` helpers) and `python/demo/pump_trip_primary_loop.py` for a
single-`PrimaryLoop` baseline-vs-perturbed comparison using `Simulation`/`postprocess.compare()`.

## References
The majority of the models were adapted from [Naghedolfeizi et al](https://trace.tennessee.edu/handle/20.500.14382/38656) and [Samet E Arda et al](https://www.researchgate.net/publication/301705924_Nonlinear_dynamic_modeling_and_simulation_of_a_passively_cooled_small_modular_reactor).


### Maintainer 
Ali Mahdi
Saad Islam
Ebny Walid Ahammed


