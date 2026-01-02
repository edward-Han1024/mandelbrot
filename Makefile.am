AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -Wall -I"$(pwd)/installroot/include"
AM_CXXFLAGS = -Wall -I"$(pwd)/installroot/include"
# libs
lib_LTLIBRARIES = lib/libmandelbrot.la
lib_libmandelbrot_la_SOURCES = mainmandelbrot.cpp mainmandelbrot.h
lib_libmandelbrot_la_LDFLAGS = -L"$(pwd)/installroot/lib" -lpng -lmpfr
