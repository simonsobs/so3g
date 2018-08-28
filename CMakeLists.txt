cmake_minimum_required (VERSION 2.8)
project (so3g)

execute_process(COMMAND mkdir -p ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
add_subdirectory(src)

