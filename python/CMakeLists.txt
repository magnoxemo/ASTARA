# =============================================================================
# nanobind Python extension module.  Built only when ASTARA_BUILD_PYTHON=ON
# (see top-level CMakeLists.txt, which locates Python and nanobind before
# add_subdirectory'ing here).  Produces `_astara*.so` in this binary
# directory; use it via `python/astara/__init__.py`, e.g.
#   PYTHONPATH=<build>/python python3 -c "import astara"
# =============================================================================

nanobind_add_module(_astara
        src/module.cpp
        src/bind_core.cpp
        src/bind_props.cpp
        src/bind_reactor.cpp
        src/bind_pump.cpp
        src/bind_pressurizer.cpp
        src/bind_sg.cpp
        src/bind_control.cpp
        src/bind_primary.cpp)
target_link_libraries(_astara PRIVATE astara::astara)

# Place the compiled extension next to the pure-Python package so
# `PYTHONPATH=<build>/python python3 -c "import astara"` works unmodified.
set_target_properties(_astara PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/astara"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/astara")

# Copy the pure-Python package (__init__.py + the network/control/sim/
# postprocess ergonomics modules) alongside the extension so the build tree
# is directly importable without an install step.
file(GLOB ASTARA_PY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/astara/*.py")
file(COPY ${ASTARA_PY_SOURCES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/astara")

# Install the compiled extension alongside the pure-Python package so
# scikit-build-core (which populates the wheel from `cmake --install`)
# includes it. DESTINATION is relative to the wheel's install prefix, so
# this lands at the "astara/" package root regardless of platform.
install(TARGETS _astara LIBRARY DESTINATION astara RUNTIME DESTINATION astara)

if(ASTARA_BUILD_TESTS)
    find_package(Python 3.8 COMPONENTS Interpreter REQUIRED)
    add_test(NAME python_smoke
             COMMAND ${Python_EXECUTABLE} -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/tests" -q)
    set_tests_properties(python_smoke PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
endif()
