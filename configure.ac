AC_INIT([mandelbrot], [1.0a1], [Edward Han])
# shellcheck disable=SC2006,SC2268 # we must support pre-POSIX shells
AM_INIT_AUTOMAKE
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([main.cpp])
: "${CXXFLAGS="-O3 -Wall"}"
AC_PROG_CXX

mkdir $srcdir/installroot

AC_LANG_PUSH([C++])

############################################################################################
#                                       ZLIB CHECKS                                        #
############################################################################################

AC_CHECK_HEADERS([zlib.h], [
    AC_CHECK_LIB([z], [deflate], [], [
        AC_MSG_NOTICE([installing zlib])
        test zlib-1.3.1.tar.gz && rm -f zlib-1.3.1.tar.gz
        test zlib-1.3.1 && rm -rf zlib-1.3.1
        wget -q "http://fresh-center.net/linux/misc/zlib-1.3.1.tar.gz" || curl -fsSLO "http://fresh-center.net/linux/misc/zlib-1.3.1.tar.gz" || exit1
        tar -xf "zlib-1.3.1.tar.gz"
        cd zlib-1.3.1 || exit 1
        AC_MSG_NOTICE([running configure])
        ./configure --prefix="$(pwd)/../installroot" > /dev/null || exit 1
        AC_MSG_NOTICE([running make])
        make -s > /dev/null || exit 1
        AC_MSG_NOTICE([running make install])
        make -s install > /dev/null || exit 1
        cd ..
    ])
], [
    AC_MSG_NOTICE([installing zlib])
    test zlib-1.3.1.tar.gz && rm -f zlib-1.3.1.tar.gz
    test zlib-1.3.1 && rm -rf zlib-1.3.1
    wget -q "http://fresh-center.net/linux/misc/zlib-1.3.1.tar.gz" || curl -fsSLO "http://fresh-center.net/linux/misc/zlib-1.3.1.tar.gz" || exit1
    tar -xf "zlib-1.3.1.tar.gz"
    cd zlib-1.3.1 || exit 1
    AC_MSG_NOTICE([running configure])
    ./configure --prefix="$(pwd)/../installroot" > /dev/null || exit 1
    AC_MSG_NOTICE([running make])
    make -s > /dev/null || exit 1
    AC_MSG_NOTICE([running make install])
    make -s install > /dev/null || exit 1
    cd ..
])

############################################################################################
#                                      LIBPNG CHECKS                                       #
############################################################################################

AC_CHECK_HEADERS([png.h], [
    AC_CHECK_LIB([png], [png_create_read_struct], [], [
        AC_MSG_NOTICE([installing libpng])
        wget -q "https://download.sourceforge.net/libpng/libpng-1.6.50.tar.gz" || curl -fsSLO "https://download.sourceforge.net/libpng/libpng-1.6.50.tar.gz" || exit 1
        tar -xf "libpng-1.6.50.tar.gz"
        cd libpng-1.6.50 || exit 1
        AC_MSG_NOTICE([running configure])
        LDFLAGS="-L../installroot/lib" CPPFLAGS="-I../installroot/include" ./configure --prefix="$(pwd)/../installroot"
        AC_MSG_NOTICE([running make])
        make -s > /dev/null || exit 1
        AC_MSG_NOTICE([running make install])
        make -s install > /dev/null || exit 1
        cd ..
    ])
], [
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
    cd ..
])

############################################################################################
#                                        GMP CHECKS                                        #
############################################################################################

AC_CHECK_HEADERS([gmp.h], [
    AC_CHECK_LIB([gmp], [__gmpz_init], [], [
        AC_MSG_NOTICE([installing gmp])
        test gmp-6.3.0.tar.xz && rm -f gmp-6.3.0.tar.xz
        test gmp-6.3.0 && rm -rf  gmp-6.3.0
        wget -q "https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz" || curl -fsSLO "https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz" || exit 1
        tar -xf "gmp-6.3.0.tar.xz"
        cd gmp-6.3.0 || exit 1
        AC_MSG_NOTICE([running configure])
        LDFLAGS="-L$(pwd)/../installroot/lib" \
        CPPFLAGS="-I$(pwd)/../installroot/include" \
        CFLAGS="-I$(pwd)/../installroot/include" \
        ./configure --prefix="$(pwd)/../installroot" > /dev/null
        AC_MSG_NOTICE([running make])
        make -s > /dev/null || exit 1
        AC_MSG_NOTICE([running make install])
        make -s install > /dev/null || exit 1
        cd ..
    ])
], [
    AC_MSG_NOTICE([installing gmp])
    test gmp-6.3.0.tar.xz && rm -f gmp-6.3.0.tar.xz
    test gmp-6.3.0 && rm -rf  gmp-6.3.0
    wget -q "https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz" || curl -fsSLO "https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz" || exit 1
    tar -xf "gmp-6.3.0.tar.xz"
    cd gmp-6.3.0 || exit 1
    AC_MSG_NOTICE([running configure])
    LDFLAGS="-L$(pwd)/../installroot/lib" \
    CPPFLAGS="-I$(pwd)/../installroot/include" \
    CFLAGS="-I$(pwd)/../installroot/include" \
    ./configure --prefix="$(pwd)/../installroot" > /dev/null
    AC_MSG_NOTICE([running make])
    make -s > /dev/null || exit 1
    AC_MSG_NOTICE([running make install])
    make -s install > /dev/null || exit 1
    cd ..
])

############################################################################################
#                                       MPFR CHECKS                                        #
############################################################################################

AC_CHECK_HEADERS([mpfr.h], [
    AC_CHECK_LIB([mpfr], [mpfr_init], [], [
        AC_MSG_NOTICE([installing mpfr])
        test mpfr-4.2.2.tar.gz && rm -f mpfr-4.2.2.tar.gz
        test mpfr-4.2.2 && rm -rf mpfr-4.2.2
        wget -q "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.gz" || curl -fsSLO "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.gz" || exit 1
        tar -xf "mpfr-4.2.2.tar.gz"
        cd mpfr-4.2.2 || exit 1
        AC_MSG_NOTICE([running configure])
        LDFLAGS="-L$(pwd)/../installroot/lib" \
        CPPFLAGS="-I$(pwd)/../installroot/include" \
        CFLAGS="-I$(pwd)/../installroot/include" \
        ./configure --prefix="$(pwd)/../installroot" > /dev/null
        AC_MSG_NOTICE([running make])
        make -s > /dev/null || exit 1
        AC_MSG_NOTICE([running make install])
        make -s install > /dev/null || exit 1
        cd ..
    ])
], [
    AC_MSG_NOTICE([installing mpfr])
    test mpfr-4.2.2.tar.gz && rm -f mpfr-4.2.2.tar.gz
    test mpfr-4.2.2 && rm -rf mpfr-4.2.2
    wget -q "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.gz" || curl -fsSLO "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.gz" || exit1
    tar -xf "mpfr-4.2.2.tar.gz"
    cd mpfr-4.2.2 || exit 1
    AC_MSG_NOTICE([running configure])
    LDFLAGS="-L$(pwd)/../installroot/lib" \
    CPPFLAGS="-I$(pwd)/../installroot/include" \
    CFLAGS="-I$(pwd)/../installroot/include" \
    ./configure --prefix="$(pwd)/../installroot" > /dev/null
    AC_MSG_NOTICE([running make])
    make -s > /dev/null || exit 1
    AC_MSG_NOTICE([running make install])
    make -s install > /dev/null || exit 1
    cd ..
])

############################################################################################
#                                       PNG++ CHECKS                                       #
############################################################################################

AC_CHECK_HEADERS([png++/png.hpp], [], [
    test png++-0.2.10.tar.gz && rm -f png++-0.2.10.tar.gz
    test png++-0.2.10 && rm -rf png++-0.2.10.tar.gz
    AC_MSG_NOTICE([installing png++ headers])
    wget -q "https://download.savannah.nongnu.org/releases/pngpp/png++-0.2.10.tar.gz" || curl -fsSLO "https://download.savannah.nongnu.org/releases/pngpp/png++-0.2.10.tar.gz" || exit 1
    tar -xf "png++-0.2.10.tar.gz"
    cd "png++-0.2.10" || exit 1
    AC_MSG_NOTICE(running make install)
    PREFIX=$(pwd)/../installroot make -s install > /dev/null || exit 1
    cd ..
])
AC_LANG_POP([C++])

mkdir $srcdir/bin
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
