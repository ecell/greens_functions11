cmake_minimum_required(VERSION 3.8)
project(gf11)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/extlib)
include_directories(${PROJECT_SOURCE_DIR}/extlib/greens_functions/build)

add_subdirectory(extlib)
add_subdirectory(src)
