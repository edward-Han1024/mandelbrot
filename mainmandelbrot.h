#include "common.h"

/**
 * @brief compute the pixel given with high precision.
 *
 * There is a long list of parameters but most of them are self-explanatory.
 * Usually there is no reason to run this function manually. It would be run
 * with mandelbrot.
 * 
 * @param temp3 units / one real pixel
 * @param temp4 units / one imaginary pixel
 */
void mandelbrot_compute_pixel (uint x, uint y, uint height, uint width,
                               uint max_iterations, uint colorlength,
                               uint precision, uint lowerprec, mpfr_t minre,
                               mpfr_t dre, mpfr_t minim, mpfr_t dim,
                               mpfr_t scale, png::image<png::rgb_pixel> &image,
                               mpfr_t one, mpfr_t two, mpfr_t bailout,
                               mpfr_t temp3, mpfr_t temp4);

/**
 * @brief Compute the pixel given with double precision.
 * There is not much reason to call this manually. This would be called by
 * mandelbrot. There is a long list of parameters but most of them are
 * self-explanatory.
 * @param temp3 units / one real pixel
 * @param temp4 units / one imaginary pixel
 */
void mandelbrot_compute_pixel (uint x, uint y, uint height, uint width,
                               uint max_iterations, uint colorlength,
                               double minre, double dre, double minim,
                               double dim, double scale,
                               png::image<png::rgb_pixel> &image,
                               double bailout, double temp3, double temp4);

/**
 * @breif Start a double precision worker.
 * There is generally no need to call this function manually. Calling mandelbrot
 * is prefered. This function has "double" in it instead of an overload as
 * std::thread doesn't support overloads.
 */
void mandelbrot_work_double (uint &x, uint &y, uint height, uint width,
                             uint max_iterations, uint colorlength,
                             double minre, double dre, double minim,
                             double dim, double scale,
                             png::image<png::rgb_pixel> &image, double bailout,
                             double temp3, double temp4);

/**
 * @breif Start a high precision worker.
 *
 * There is generally no need to call this manually. Calling mandelbrot is
 * prefered.
 */
void mandelbrot_work (uint &x, uint &y, uint height, uint width,
                      uint max_iterations, uint colorlength, uint precision,
                      uint lowerprec, mpfr_t minre, mpfr_t dre, mpfr_t minim,
                      mpfr_t dim, mpfr_t scale,
                      png::image<png::rgb_pixel> &image, mpfr_t one,
                      mpfr_t two, mpfr_t bailout, mpfr_t temp3, mpfr_t temp4);

/**
 * @breif Generate the mandelbrot fractal given the parameters using arbirtrary
 * precision arithmetic.
 *
 * This function currently uses the mpfr library to do arbitrary precision
 * arithmetic. If the necessary precision is 53 bits or lower, please use the
 * other function instead, as that will be much faster.
 
 * This specific function currently does not use bulb checking for the cardiod
 * and the bulb to the left. It currently uses the basic escape speed method.
 * The coloring method achieved by this function is a hsv gradient with h adjusted,
 * and s and v constant. The colorlength is how many iterations it would take to
 * fully cycle through the palette.
 * 
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param precision (u_long) the precision to use for the calculation
 * @param lowprec (u_long) the precision for things that don't need to be as
 * precise such as the aspect ration or scale. Generally, low values for this
 * could also make the image blocky
 * @param lowerprec (u_long) the precision for things that are integers but
 * need a mpfr_t value so we can access it easily. Should be at least
 * log_2(width)
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
int mandelbrot (uint height, uint width, u_long precision, u_long loprec,
		u_long lowerprec, const char *CENTERRE, const char *CENTERIM,
		const char *BAILOUT, const char *SCALE, uint max_iterations,
		uint colorlength);

/**
 * @brief Generate the mandelbrot fractal given parameters, but without
 * arbitrary precision support.
 * 
 * This function uses double precision, which is equivalent to 53 bits in mpfr.
 * Thus, if the precision is set to less than 53 in the arbitrary precision
 * function, then this overload should be used instead.
 *
 * This specific function currently uses bulb checking for the cardiod and the
 * bulb to the left. It currently uses the basic escape speed method. The
 * coloring method achieved by this function is a hsv gradient with h adjusted, and
 * s and v constant. The colorlength is how many iterations it would take to fully
 * cycle through the palette.
 *
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param centerre (double) the center of the real axis
 * @param centerim (double) the center of the imaginary axis
 * @param bailout (double)  the square of the distance away from the origin needed to be considered as diverged
 * @param scale (double) the scale of the image (zoom)
 * @param max_iterations the maximum number of iterations
 * @param colorlength the length of the palette
 * @returns (int) the status (exit value) 0 is no errors
 */
int mandelbrot (uint height, uint width, double centerre, double centerim,
                double bailout, double scale, uint max_iterations,
                uint colorlength);

