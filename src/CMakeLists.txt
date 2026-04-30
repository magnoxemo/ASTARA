# =============================================================================
# Per-module static libraries.  Each is small and self-contained so individual
# components can be linked into focused tests.
# =============================================================================

# ---- core: header-only (Units, State, Integrator, Port) ----------------------
add_library(astara_core INTERFACE)
target_include_directories(astara_core INTERFACE
        "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_core INTERFACE astara_warnings)
if(ASTARA_ENABLE_OPENMP)
    target_link_libraries(astara_core INTERFACE OpenMP::OpenMP_CXX)
endif()

# ---- props: WaterProperties + IF97Water + LinearizedWater --------------------
add_library(astara_props
        props/IF97Water.cpp
        props/LinearizedWater.cpp)
target_include_directories(astara_props PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_props
        PUBLIC  astara_core
        PRIVATE astara_if97
                astara_warnings)

# ---- reactor: PointKinetics + ReactorThermal + Reactor -----------------------
add_library(astara_reactor
        reactor/PointKinetics.cpp
        reactor/ReactorThermal.cpp
        reactor/Reactor.cpp)
target_include_directories(astara_reactor PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_reactor
        PUBLIC  astara_core
        PRIVATE astara_warnings)

# ---- pump: HomologousPump ----------------------------------------------------
add_library(astara_pump
        pump/HomologousPump.cpp)
target_include_directories(astara_pump PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_pump
        PUBLIC  astara_core
        PRIVATE astara_warnings)

# ---- pressurizer: Pressurizer ------------------------------------------------
add_library(astara_pressurizer
        pressurizer/Pressurizer.cpp)
target_include_directories(astara_pressurizer PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_pressurizer
        PUBLIC  astara_core
                astara_props
        PRIVATE astara_warnings)

# ---- sg: Ali + nodal steam generators ----------------------------------------
add_library(astara_sg
        sg/AliSteamGenerator.cpp
        sg/HelicalCoilSteamGenerator.cpp)
target_include_directories(astara_sg PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_sg
        PUBLIC  astara_core
                astara_props
        PRIVATE astara_warnings)

# ---- control: PID + reactor + pressurizer + 3-element ------------------------
# Header-only.
add_library(astara_control INTERFACE)
target_include_directories(astara_control INTERFACE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_control INTERFACE astara_warnings)

# ---- primary: integrated loop -----------------------------------------------
add_library(astara_primary
        primary/PrimaryLoop.cpp)
target_include_directories(astara_primary PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(astara_primary
        PUBLIC  astara_core
                astara_reactor
                astara_sg
                astara_pump
                astara_pressurizer
                astara_control
        PRIVATE astara_warnings)
if(ASTARA_ENABLE_OPENMP)
    target_link_libraries(astara_primary PRIVATE OpenMP::OpenMP_CXX)
endif()

# =============================================================================
# Umbrella interface library.  Most users -- and most of the tests / examples
# -- should link to this rather than picking individual sub-libraries.
# =============================================================================
add_library(astara INTERFACE)
target_link_libraries(astara INTERFACE
        astara_core
        astara_props
        astara_reactor
        astara_pump
        astara_pressurizer
        astara_sg
        astara_control
        astara_primary)
add_library(astara::astara ALIAS astara)
