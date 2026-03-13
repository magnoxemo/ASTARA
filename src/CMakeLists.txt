file(GLOB src *.cpp *.C)

add_library(astara_lib SHARED ${src})
target_link_libraries(astara_lib
        MPI::MPI_CXX
        CoolProp
)

target_include_directories(astara_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)