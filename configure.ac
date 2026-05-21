AC_INIT([mandelbrot], [1.0a1], [Edward Han])
# shellcheck disable=SC2006,SC2268 # we must support pre-POSIX shells
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([mainmandelbrot.cpp])
: "${CXXFLAGS="-O3 -Wall"}"
AC_PROG_CXX
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

test lib || mkdir lib
test bin || mkdir bin
test local || mkdir local
