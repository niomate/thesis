cmake_minimum_required(VERSION 3.13)
project(thesis)

set(CMAKE_C_STANDARD 11)
add_definitions("-Wall" "-g" "-DCMAKE_BUILD_TYPE=Debug")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin/)
# Link math library
link_libraries(m)

file(GLOB AMSS_SRC src/amss/*.c)

add_executable(test src/amss/amss_corner_detection.c src/amss/amss.c src/amss/chain.c src/amss/least_squares.c src/utils.c src/test_images.c)
add_executable(amss_corner ${AMSS_SRC} src/utils.c)
add_executable(inpainting src/inpainting/inpainting.c)
add_executable(corner_detection src/corner_detection.c src/utils.c src/corner_detection_main.c)

add_library(amss  SHARED ${AMSS_SRC} src/utils.c)
add_library(inpaint SHARED src/inpainting/inpainting.c)
add_library(mask SHARED src/corner_detection.c src/mask_algo.c src/utils.c)

