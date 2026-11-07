cmake_minimum_required(VERSION 3.16)
include(CheckIncludeFiles)

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
check_include_files(gmp.h GMP_H)
if(NOT GMP_H)
message(WARNING GMP is a required dependency of MPFR. The header gmp.h was not found.)
message(WARNING This most likely means GMP and MPFR wasn't installed correctly.)
endif()
check_include_files(mpfr.h MPFR_H)
if(NOT MPFR_H)
message(SEND_ERROR MPFR is a required dependency. The header mpfr.h was not found.)
endif()
check_include_files(boost/math/interpolators/pchip.hpp BOOST_PCHIP_HPP)
if(NOT BOOST_PCHIP_HPP)
message(SEND_ERROR Boost is a required dependency. The header boost/math/interpolators/pchip.hpp was not found.)
endif()
check_include_files(png++/png.hpp PNG_HPP)
if(NOT PNG_HPP)
message(SEND_ERROR PNG++ is a required dependency. The header png++/png.hpp was not found.)
endif()


# targets
add_library(mandelbrot SHARED
    src/mandelbrot.cpp
    src/util.cpp
    src/burningship.cpp
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