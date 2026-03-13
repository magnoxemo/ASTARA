cmake_minimum_required(VERSION 3.11)
project(astara)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


find_package(MPI REQUIRED)
find_package(GTest REQUIRED)

set(COOLPROP_STATIC_LIBRARY ON CACHE BOOL "Build CoolProp as static library")
set(COOLPROP_SHARED_LIBRARY OFF CACHE BOOL "Don't build shared library")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/externals/CoolProp")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/externals/CoolProp/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/externals/CoolProp/externals/fmtlib)

file(GLOB src ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB unit_test_src ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/*.cpp)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/unit_tests")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/examples")