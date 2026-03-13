cmake_minimum_required(VERSION 3.11)
project(astara)

find_package(MPI REQUIRED)
find_package(GTest REQUIRED)

set(CMAKE_BUILD_TYPE Release)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")


file(GLOB src ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB unit_test_src ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/*.cpp)


add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/unit_tests")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/examples")