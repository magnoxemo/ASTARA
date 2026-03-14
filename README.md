# ASTARA

**ASTARA** is a simulation tool for modeling the **primary coolant circuit dynamics of a Pressurized Water Reactor (PWR)** nuclear power plant.
The project began as a personal hobby during my undergraduate studies and has since evolved into an ongoing development effort. 
The current implementation is written in **C++**, and the internal API is still evolving. As a result, **breaking API changes may occur frequently**.

ASTARA uses [**CoolProp**](https://github.com/CoolProp/CoolProp) for thermodynamic property lookups.


## Features

- Simulation of **PWR primary coolant loop dynamics**
- OOP **C++ architecture**
- MPI for parallelism (that's a todo) 
- Depends on **CoolProp** for coolant thermal property lookup (water)
- **CMake-based** build system

## Installation

### Just copy the commands and run one by one

```bash
git clone https://github.com/magnoxemo/ASTARA.git
cd ASTARA

cd externals
git submodule update --init --recursive
cd ..

mkdir build
cd build
cmake ..
make -j $(nproc)
```
