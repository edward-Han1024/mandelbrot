# mandelbrot

Mandelbrot is a library to render some specific fractals.

It currently uses `cmake` to build. It depends on GNU MPFR, PNG++, and boost.
GNU MPFR and its dependency GNU GMP is licensed under the GNU GPL, and is not included. (you must install it)
PNG++ and it's dependencies ZLIB and libpng are under the 3 clause BSD license, the ZLIB license, and the libpng license, respectively. They are also not included.
Boost is under the Boost license, and is not included. The specific header needed is `boost/math/interpolators/pchip.hpp` for gradient interpolation.