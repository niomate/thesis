cmake_minimum_required(VERSION 3.13)
project(thesis)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_C_STANDARD 11)

# Link math library
link_libraries(m)

add_executable(corners src/corner_detection.c src/utils.c)
add_executable(inpainting src/inpainting.c)

add_library(corner SHARED src/corner_detection.c src/utils.c)
add_library(inpaint SHARED src/inpainting.c)

