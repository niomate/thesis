cmake_minimum_required(VERSION 3.13)
project(thesis)

set(CMAKE_C_STANDARD 11)
add_definitions("-Wall" "-g" "-DCMAKE_BUILD_TYPE=Debug")

# Link math library
link_libraries(m)

file(GLOB AMSS_SRC src/amss/*.c)
#file(GLOB STRUCT_SRC src/structure_tensor/*.c)

add_executable(test src/amss/amss_corner_detection.c src/amss/amss.c src/amss/chain.c src/amss/least_squares.c src/utils.c src/test/test_images.c)
add_executable(amss_corner ${AMSS_SRC} src/utils.c)
add_executable(inpainting src/inpainting/inpainting.c)
#add_executable(corner_detection ${STRUCT_SRC} src/utils.c)
#
#
#add_library(amss  SHARED ${AMSS_SRC} src/utils.c)
#add_library(inpaint SHARED src/inpainting/inpainting.c)
#add_library(mask SHARED src/structure_tensor/corner_detection.c src/mask/main.c src/utils.c)
#
