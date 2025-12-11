include_directories(.)
include_directories(../common)
include_directories(../src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin) 

add_executable(test_ctru_compact test_kem.c ../src/reduce.c)
target_link_libraries(test_ctru_compact src)