AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4
# libs
lib_LTLIBRARIES = lib/libmandelbrot.la
lib_libmandelbrot_la_SOURCES = mainmandelbrot.cpp util.cpp mainburningship.cpp
lib_libmandelbrot_la_LDFLAGS = -lgmp -lpng -lmpfr

bin_PROGRAMS = bin/main
bin_main_SOURCES = main.cpp
bin_main_LDADD = lib/libmandelbrot.la
bin_main_LDFLAGS = -lgmp -lpng -lmpfr
