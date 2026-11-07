// mandelbrot.h: the file which has the public api
// Copyright (C) 2025  Edward Han
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once
#ifndef MANDELBROT_H
#define MANDELBROT_H
#include <cstdint>
#include <string>

#define uint uint64_t
#define u_long uint
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

/**
 * @brief Generate the burning ship fractal given the parameters using
 * arbirtrary precision arithmetic.
 *
 * This function currently uses the mpfr library to do arbitrary precision
 * arithmetic. If the necessary precision is 53 bits or lower, please use the
 * other function instead, as that will be much faster.
 *
 * This specific function currently does not use bulb checking for the cardiod
 * and the bulb to the left. It currently uses the basic escape speed method.
 * The coloring method achieved by this function is a hsv gradient with h
 adjusted,
 * and s and v constant. The colorlength is how many iterations it would take
 * to fully cycle through the palette.
 *
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param precision (u_long) the precision to use for the calculation
 * @param lowprec (u_long) the precision for things that don't need to be as
 * precise such as the aspect ration or scale. Generally, low values for this
 * could also make the image blocky
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
int burningship(uint64_t height, uint64_t width, u_long precision,
                u_long loprec, u_long lowerprec, const char *CENTERRE,
                const char *CENTERIM, const char *BAILOUT, const char *SCALE,
                uint64_t max_iterations, uint64_t colorlength);

/**
 * @brief Generate the burning ship fractal given parameters, but without
 * arbitrary precision support.
 *
 * This function uses double precision, which is equivalent to 53 bits in mpfr.
 * Thus, if the precision is set to less than 53 in the arbitrary precision
 * function, then this overload should be used instead.
 *
 * This specific function currently uses bulb checking for the cardiod and the
 * bulb to the left. It currently uses the basic escape speed method. The
 * coloring method achieved by this function is a hsv gradient with h adjusted,
 * and s and v constant. The colorlength is how many iterations it would take
 * to fully cycle through the palette.
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
int burningship(uint64_t height, uint64_t width, double centerre,
                double centerim, double bailout, double scale,
                uint64_t max_iterations, uint64_t colorlength,
                std::string filename);
#endif
