file(GLOB unit_test_src *.cpp *.C)

#adding test executables for every case
foreach(TEST_SOURCE ${unit_test_src})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} astara_lib MPI::MPI_CXX GTest::GTest GTest::Main )

endforeach()