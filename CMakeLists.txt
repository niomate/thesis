cmake_minimum_required(VERSION 3.13)
project(thesis)

set(CMAKE_C_STANDARD 11)
add_definitions("-Wall" "-g" "-DCMAKE_BUILD_TYPE=Debug")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin/)
# Link math library
link_libraries(m)

file(GLOB AMSS_SRC src/amss/*.c)

add_executable(amss_corner ${AMSS_SRC} src/utils.c)
add_executable(inpainting src/inpainting/inpainting.c)
add_executable(corner_detection src/corner_detection.c src/utils.c src/corner_detection_main.c)

add_library(libamss.so  SHARED ${AMSS_SRC})
add_library(libinpaint.so SHARED src/inpainting/inpainting.c)
add_library(libmask.so SHARED src/corner_detection.c src/mask_algo.c src/utils.c)

