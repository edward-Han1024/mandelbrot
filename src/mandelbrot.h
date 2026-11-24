#include "common.h"
#pragma once

/**
 * @brief Generate the mandelbrot fractal given the parameters using arbitrary
 * precision arithmetic.
 *
 * This function currently uses the mpfr library to do arbitrary precision
 * arithmetic. If the necessary precision is 53 bits or lower, the other
 * function is faster.
 *
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param precision (u_long) the precision to use for the calculation
 * @param CENTERRE (const char *) the c-style string for the center of the real
 * axis
 * @param CENTERIM (const char *) the c-style string for the center of the
 * imaginary axis
 * @param BAILOUT (const char *) the c-style string for the square of the
 * distance away from the origin needed to be considered as diverged
 * @param SCALE (const char *) the c-style string for the scale of the image
 * (zoom)
 * @param max_iterations the maximum number of iterations
 * @param colorlength the length of the palette
 * @returns (int) the status (exit value) 0 is no errors
 */

int mandelbrot(uint64_t height, uint64_t width, u_long precision,
               const char *CENTERRE, const char *CENTERIM, const char *BAILOUT,
               const char *SCALE, uint64_t max_iterations, uint64_t colorlength,
               std::string filename);

/**
 * @brief Generate the mandelbrot fractal given parameters, but without
 * arbitrary precision support.
 *
 * This function uses double precision, which is equivalent to 53 bits in
 * mpfr. Thus, if the precision is set to less than 53 in the arbitrary
 * precision function, then this overload should be used instead.
 *
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param centerre (double) the center of the real axis
 * @param centerim (double) the center of the imaginary axis
 * @param bailout (double)  the square of the distance away from the origin
 * needed to be considered as diverged
 * @param scale (double) the scale of the image (zoom)
 * @param max_iterations the maximum number of iterations
 * @param colorlength the length of the palette
 * @returns (int) the status (exit value) 0 is no errors
 */
int mandelbrot(uint64_t height, uint64_t width, double centerre,
               double centerim, double bailout, double scale,
               uint64_t max_iterations, uint64_t colorlength,
               std::string filename);
