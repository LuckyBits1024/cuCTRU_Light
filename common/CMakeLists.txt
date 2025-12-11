# 添加头文件的搜索路径
include_directories(.)

aux_source_directory(. COMMON_LIST)

add_library(common STATIC ${COMMON_LIST}) # 生成静态库
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib) # 将静态库放在指定目录
