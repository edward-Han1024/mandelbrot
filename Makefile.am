AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -Wall -I"$(shell pwd)/gmp-6.3.0 -I$(shell pwd)/installroot/include"
AM_CXXFLAGS = -Wall -I"$(shell pwd)/gmp-6.3.0 -I$(shell pwd)/installroot/include"
AM_CFLAGS = -Wall -I"$(shell pwd)/gmp-6.3.0 -I$(shell pwd)/installroot/include"
AM_LDFLAGS = -Wall -L"$(shell pwd)/gmp-6.3.0 -L$(shell pwd)/installroot/include"
SUBDIRS = gmp-6.3.0 mpfr-4.2.2 .
# libs
lib_LTLIBRARIES = lib/libmandelbrot.la
lib_libmandelbrot_la_SOURCES = mainmandelbrot.cpp mainmandelbrot.h
lib_libmandelbrot_la_LDFLAGS = -lpng -lmpfr
bin_PROGRAMS = bin/main
bin_main_SOURCES = main.cpp
bin_main_LDFLAGS = -Wall -L"$(shell pwd)/lib" -lmandelbrot
