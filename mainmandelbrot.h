#include "common.h"

/**
 * @brief Convertes hsv value into an rgb value.
 * @param h (double) The hue value. The range is [0, 1), but it can loop
 * around.
 * @param s (double) The saturation of the color. The range is [0, 1)
 * @param v (double) The value of the color. The range is [0, 1)
 * @param r (uint&) A pointer to the uint to store the value of r. The range will be [0, 255].
 * @param g (uint&) A pointer to the uint to store the value of g. The range will be [0, 255].
 * @param b (uint&) A pointer to the uint to store the value of b. The range
 * will be [0, 255].
 */
void hsvToRgb (double h, double s, double v, uint &r, uint &g, uint &b);

/**
 * @brief compute the pixel given with high precision.
 * There is a long list of parameters but most of them are self-explanatory. Usually there is no reason to run this function manually. It would be run with maindelbrot.
 * \c temp3 = units / one real pixel
 * \c temp4 = units / one imaginary pixel
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
 * There is not much reason to call this manually. This would be called by mandelbrot. There is a long list of parameters but most of them are self-explanatory.
 * \c temp3 = units / one real pixel
 * \c temp4 = units / one imaginary pixel
 */
void mandelbrot_compute_pixel (uint x, uint y, uint height, uint width,
                               uint max_iterations, uint colorlength,
                               double minre, double dre, double minim,
                               double dim, double scale,
                               png::image<png::rgb_pixel> &image,
                               double bailout, double temp3, double temp4);

/**
 * @deprecated
 * This function will be removed soon and a replacement made using thread pools.
 */
void burningship_compute_pixels (
				 uint start_x, uint end_x, uint height, uint width, uint max_iterations,
				 uint colorlength, uint precision, uint lowerprec, mpfr_t minre, mpfr_t dre,
				 mpfr_t minim, mpfr_t dim, mpfr_t scale, png::image<png::rgb_pixel> &image,
				 mpfr_t one, mpfr_t two, mpfr_t bailout, mpfr_t temp3, mpfr_t temp4);


/**
 * the function to generate the mandelbrot fractal given the parameters
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
 * @brief the function to generate the mandelbrot fractal given the
 * parameters, but without arbitrary precision support
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param centerre (double) the center of the real axis
 * @param centerim (double the center of the imaginary axis
 * @param bailout (double)  the square of the distance away from the
 * origin needed to be considered as diverged
 * @param scale (double) the scale of the image (zoom)
 * @param max_iterations the maximum number of iterations
 * @param colorlength the length of the palette
 * @returns (int) the status (exit value) 0 is no errors
 */
int mandelbrot (uint height, uint width, double centerre, double centerim,
                double bailout, double scale, uint max_iterations,
                uint colorlength);

/**
 * @brief the function to generate the mandelbrot fractal given the parameters
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
int burningship (uint height, uint width, u_long precision, u_long loprec,
                 u_long lowerprec, const char *CENTERRE, const char *CENTERIM,
                 const char *BAILOUT, const char *SCALE, uint max_iterations,
                 uint colorlength);

void mandelbrot_work_double (uint &x, uint &y, uint height, uint width,
                             uint max_iterations, uint colorlength,
                             double minre, double dre, double minim,
                             double dim, double scale,
                             png::image<png::rgb_pixel> &image, double bailout,
                             double temp3, double temp4);

void mandelbrot_work (uint &x, uint &y, uint height, uint width,
                      uint max_iterations, uint colorlength, uint precision,
                      uint lowerprec, mpfr_t minre, mpfr_t dre, mpfr_t minim,
                      mpfr_t dim, mpfr_t scale,
                      png::image<png::rgb_pixel> &image, mpfr_t one,
                      mpfr_t two, mpfr_t bailout, mpfr_t temp3, mpfr_t temp4);
