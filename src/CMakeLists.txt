file(GLOB src *.cpp *.C)

add_library(astara_lib SHARED ${src})
target_link_libraries( astara_lib MPI::MPI_CXX)