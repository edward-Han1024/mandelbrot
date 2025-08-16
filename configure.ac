AC_INIT([mandelbrot], [1.0a1], [Edward Han])
# shellcheck disable=SC2006,SC2268 # we must support pre-POSIX shells
# also this is ./missing thing but here it is for less notes
AC_CONFIG_SRCDIR([main.cpp])
: "${CXXFLAGS="-O3 -Wall"}"
AC_PROG_CXX

mkdir $srcdir/installroot
# libpng checks for a function
AC_CHECK_LIB([png], [png_create_read_struct], [], [
    AC_MSG_NOTICE([installing libpng])
    wget -q "https://download.sourceforge.net/libpng/libpng-1.6.50.tar.gz" || curl -fsSLO "https://download.sourceforge.net/libpng/libpng-1.6.50.tar.gz" || exit1
    tar -xf "libpng-1.6.50.tar.gz"
    cd libpng-1.6.50 || exit 1
    AC_MSG_NOTICE([running configure])
    ./configure --prefix="$(pwd)/../installroot" > /dev/null || exit 1
    AC_MSG_NOTICE([running make])
    make -s > /dev/null || exit 1
    AC_MSG_NOTICE([running make install])
    make -s install > /dev/null || exit 1
])

# mpfr checks for a function
AC_CHECK_LIB([mpfr], [mpfr_init], [], [
    AC_MSG_NOTICE([installing mpfr])
    wget -q "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.gz" || curl -fsSLO "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.gz" || exit1
    tar -xf "mpfr-4.2.2.tar.gz"
    cd mpfr-4.2.2 || exit 1
    AC_MSG_NOTICE([running configure])
    ./configure --prefix="$(pwd)/../installroot" > /dev/null || exit 1
    AC_MSG_NOTICE([running make])
    make -s > /dev/null || exit 1
    AC_MSG_NOTICE([running make install])
    make -s install > /dev/null || exit 1
])

# png++ checks for a HEADER
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([png++/png.hpp], [], [
    AC_MSG_NOTICE([installing png++ headers])
    wget -q "https://download.savannah.nongnu.org/releases/pngpp/png++-0.2.10.tar.gz" || curl -fsSLO "https://download.savannah.nongnu.org/releases/pngpp/png++-0.2.10.tar.gz" || exit 1
    tar -xf "png++-0.2.10.tar.gz"
    cd "png++-0.2.10" || exit 1
    AC_MSG_NOTICE(running make install)
    PREFIX=$(pwd)/../installroot make -s install > /dev/null || exit 1
])
AC_LANG_POP([C++])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
