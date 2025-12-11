include_directories(.)
include_directories(../common)

aux_source_directory(. SRC_LIST)

add_library(poly STATIC poly.c ntt.c inverse.c)
add_library(src STATIC ${SRC_LIST})
target_link_libraries(src common)
