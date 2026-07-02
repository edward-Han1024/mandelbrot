cmake_minimum_required(VERSION 3.16)

project(mandelbrot VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
endif()

find_package(PNG REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(MPFR_LIBRARY mpfr REQUIRED)

# targets
add_library(mandelbrot SHARED
    src/mainmandelbrot.cpp
    src/util.cpp
    src/mainburningship.cpp
)

target_sources(mandelbrot
    PUBLIC FILE_SET headers
    TYPE HEADERS
    BASE_DIRS "include"
    FILES 
        include/mandelbrot.h
)

target_link_libraries(mandelbrot
    PRIVATE
        ${GMP_LIBRARY}
        PNG::PNG
        ${MPFR_LIBRARY}
)

add_executable(main
    src/main.cpp
)

target_link_libraries(main
    PRIVATE
        mandelbrot
)

install(TARGETS mandelbrot
    LIBRARY DESTINATION lib
    FILE_SET headers DESTINATION include
)

install(TARGETS main
    RUNTIME DESTINATION bin
)