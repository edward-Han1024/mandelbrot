#include"common.h"

/**
 * @brief Convertes hsv value into an rgb value.
 * @param h (double) The hue value. The range is [0, 1), but it can loop around.
 * @param s (double) The saturation of the color. The range is [0, 1)
 * @param v (double) The value of the color. The range is [0, 1)
 * @param r (uint&) A pointer to the uint to store the value of r. The range will be [0, 255].
 * @param g (uint&) A pointer to the uint to store the value of g. The range will be [0, 255].
 * @param b (uint&) A pointer to the uint to store the value of b. The range will be [0, 255].
 */
void hsvToRgb(double h, double s, double v, uint &r, uint &g, uint &b);

/**
 * @brief Computes the set of pixels given.
 * There is a long list of parameters but most of them are self-explanatory
 * \c temp3 = units / one real pixel
 * \c temp4 = units / one imaginary pixel
 */
void mandelbrot_compute_pixels(uint start_x, uint end_x, uint height, uint width,
        uint max_iterations, uint colorlength, uint precision, uint lowerprec,
        mpfr_t minre, mpfr_t dre, mpfr_t minim, mpfr_t dim, 
        mpfr_t scale, png::image<png::rgb_pixel> &image,
        mpfr_t one, mpfr_t two, mpfr_t bailout, mpfr_t temp3, mpfr_t temp4
);

/**
 * @brief Computes the set of pixels given.
 * There is a long list of parameters but most of them are self-explanatory
 * \c temp3 = units / one real pixel
 * \c temp4 = units / one imaginary pixel
 */
void burningship_compute_pixels(uint start_x, uint end_x, uint height, uint width,
        uint max_iterations, uint colorlength, uint precision, uint lowerprec,
        mpfr_t minre, mpfr_t dre, mpfr_t minim, mpfr_t dim, 
        mpfr_t scale, png::image<png::rgb_pixel> &image,
        mpfr_t one, mpfr_t two, mpfr_t bailout, mpfr_t temp3, mpfr_t temp4
);

/**
 * @brief the function to generate the mandelbrot fractal given the parameters
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param precision (u_long) the precision to use for the calculation
 * @param lowprec (u_long) the precision for things that don't need to be as precise such as the aspect ration or scale. Generally, low values for this could also make the image blocky
 * @param lowerprec (u_long) the precision for things that are integers but need a mpfr_t value so we can access it easily. Should be at least log_2(width)
 * @param CENTERRE (const char *) the c-style string for the center of the real axis
 * @param CENTERIM (const char *) the c-style string for the center of the imaginary axis
 * @param BAILOUT (const char *) the c-style string for the square of the distance away from the origin needed to be considered as diverged
 * @param SCALE (const char *) the c-style string for the scale of the image (zoom)
 * @param max_iterations the maximum number of iterations
 * @param colorlength the length of the palette
 * @returns (int) the status (exit value) 0 is no errors
 */
int mandelbrot(uint height, uint width,
    u_long precision, u_long loprec, u_long lowerprec,
    const char * CENTERRE, const char * CENTERIM,
    const char * BAILOUT, const char * SCALE,
    uint max_iterations, uint colorlength
);

/**
 * @brief the function to generate the mandelbrot fractal given the parameters
 * @param height (uint) the height of the image
 * @param width (uint) the width of the image
 * @param precision (u_long) the precision to use for the calculation
 * @param lowprec (u_long) the precision for things that don't need to be as precise such as the aspect ration or scale. Generally, low values for this could also make the image blocky
 * @param lowerprec (u_long) the precision for things that are integers but need a mpfr_t value so we can access it easily. Should be at least log_2(width)
 * @param CENTERRE (const char *) the c-style string for the center of the real axis
 * @param CENTERIM (const char *) the c-style string for the center of the imaginary axis
 * @param BAILOUT (const char *) the c-style string for the square of the distance away from the origin needed to be considered as diverged
 * @param SCALE (const char *) the c-style string for the scale of the image (zoom)
 * @param max_iterations the maximum number of iterations
 * @param colorlength the length of the palette
 * @returns (int) the status (exit value) 0 is no errors
 */
int burningship(uint height, uint width,
    u_long precision, u_long loprec, u_long lowerprec,
    const char * CENTERRE, const char * CENTERIM,
    const char * BAILOUT, const char * SCALE,
    uint max_iterations, uint colorlength
);
