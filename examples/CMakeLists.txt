file(GLOB examples *.cpp *.C)

#adding test executables for every case
foreach(EXAMPLE_SOURCE ${examples})
    get_filename_component(example_name ${EXAMPLE_SOURCE} NAME_WE)
    add_executable(${example_name} ${EXAMPLE_SOURCE})
    target_link_libraries(${example_name} astara_lib MPI::MPI_CXX GTest::GTest GTest::Main )

endforeach()