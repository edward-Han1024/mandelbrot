// mainmandelbrot: the file describing the mandelbrot and burningship fractals
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

#include"common.h"

std::mutex image_mutex;

void hsvToRgb(double h, double s, double v, uint &r, uint &g, uint &b) {
    /**
     * @brief Convertes hsv value into an rgb value.
     * @param h (double) The hue value. The range is [0, 1), but it can loop around.
     * @param s (double) The saturation of the color. The range is [0, 1)
     * @param v (double) The value of the color. The range is [0, 1)
     * @param r (uint&) A pointer to the uint to store the value of r. The range will be [0, 255].
     * @param g (uint&) A pointer to the uint to store the value of g. The range will be [0, 255].
     * @param b (uint&) A pointer to the uint to store the value of b. The range will be [0, 255].
     */
    double c = v * s;
    double hPrime = fmod(6.0 * h, 6.0); // Account for looping over
    double x = c * (1 - std::abs(fmod(hPrime, 2) - 1));
    double r1, g1, b1;
    // There are 6 cases for the conversion
    if (hPrime >= 0 && hPrime < 1) {
        r1 = c; g1 = x; b1 = 0;
    } else if (hPrime >= 1 && hPrime < 2) {
        r1 = x; g1 = c; b1 = 0;
    } else if (hPrime >= 2 && hPrime < 3) {
        r1 = 0; g1 = c; b1 = x;
    } else if (hPrime >= 3 && hPrime < 4) {
        r1 = 0; g1 = x; b1 = c;
    } else if (hPrime >= 4 && hPrime < 5) {
        r1 = x; g1 = 0; b1 = c;
    } else { // hPrime >= 5 && hPrime < 6
        r1 = c; g1 = 0; b1 = x;
    }
    // brightness correction
    double m = v - c;
    r = (r1 + m) * 255.0;
    g = (g1 + m) * 255.0;
    b = (b1 + m) * 255.0;
}

void mandelbrot_compute_pixels(uint start_x, uint end_x, uint height, uint width,
        uint max_iterations, uint colorlength, uint precision, uint lowerprec,
        mpfr_t minre, mpfr_t dre, mpfr_t minim, mpfr_t dim, 
        mpfr_t scale, png::image<png::rgb_pixel> &image,
        mpfr_t one, mpfr_t two, mpfr_t bailout, mpfr_t temp3, mpfr_t temp4)
{
    /**
     * @brief Computes the set of pixels given.
     * There is a long list of parameters but most of them are self-explanatory
     * \c temp3 = units / one real pixel
     * \c temp4 = units / one imaginary pixel
     */
    // init r g and b
    uint r = 0, g = 0, b = 0;
    // initialize variables
    mpfr_t px, py, temp, temp2, re, im, wre, wim;
    mpfr_init2(px, lowerprec); mpfr_init2(py, lowerprec); 
    mpfr_init2(temp, precision); mpfr_init2(temp2, precision);
    mpfr_init2(re, precision); mpfr_init2(im, precision); 
    mpfr_init2(wre, precision); mpfr_init2(wim, precision);
    // this is the loop itself
    for (uint x = start_x; x < end_x; ++x) {
        mpfr_set_d(px, x, MPFR_RNDN); // a mpfr_t x value
        for (uint y = 0; y < height; ++y) {
            mpfr_set_d(py, y, MPFR_RNDN); // a mpfr_t y value
            mpfr_mul(temp, temp3, px, MPFR_RNDN); mpfr_add(re, minre, temp, MPFR_RNDN); // calculate re
            mpfr_mul(temp, temp4, py, MPFR_RNDN); mpfr_add(im, minim, temp, MPFR_RNDN); // calculate im
            mpfr_set_d(wre, 0, MPFR_RNDN); mpfr_set_d(wim, 0, MPFR_RNDN);
            uint t = 0; // the time
            while (t < max_iterations) {
                // do the real calc
                mpfr_sqr(temp, wre, MPFR_RNDN); mpfr_sqr(temp2, wim, MPFR_RNDN);
                mpfr_sub(temp, temp, temp2, MPFR_RNDN);
                mpfr_add(temp, temp, re, MPFR_RNDN); // add term
                // do the imaginary calc
                mpfr_mul(temp2, wim, wre, MPFR_RNDN); mpfr_mul(temp2, temp2, two, MPFR_RNDN);
                mpfr_add(temp2, temp2, im, MPFR_RNDN); // add term
                mpfr_set(wre, temp, MPFR_RNDN); mpfr_set(wim, temp2, MPFR_RNDN); // set values
                // calculate bailout
                mpfr_mul(temp, wre, wre, MPFR_RNDN); mpfr_mul(temp2, wim, wim, MPFR_RNDN);
                mpfr_add(temp, temp, temp2, MPFR_RNDN);
                if (mpfr_cmp(temp, bailout) > 0) break;
                ++t;
            }
            // store pixels
            if (t == max_iterations) {
                std::lock_guard<std::mutex> lock(image_mutex);
                image[height - y - 1][x] = png::rgb_pixel(0, 0, 0);
                continue;
            }
            hsvToRgb(((double) t) / colorlength, 1, 1, r, g, b);
            {
                std::lock_guard<std::mutex> lock(image_mutex);
                image[height - y - 1][x] = png::rgb_pixel(r, g, b);
            }
        }
    }
    // clear the (small amount of) mpfr_t variables we defined
    mpfr_clear(px); mpfr_clear(py); mpfr_clear(temp); mpfr_clear(temp2);
    mpfr_clear(re); mpfr_clear(im);
    mpfr_clear(wre); mpfr_clear(wim);
}

void burningship_compute_pixels(uint start_x, uint end_x, uint height, uint width,
        uint max_iterations, uint colorlength, uint precision, uint lowerprec,
        mpfr_t minre, mpfr_t dre, mpfr_t minim, mpfr_t dim, 
        mpfr_t scale, png::image<png::rgb_pixel> &image,
        mpfr_t one, mpfr_t two, mpfr_t bailout, mpfr_t temp3, mpfr_t temp4)
{
    /**
     * @brief Computes the set of pixels given.
     * There is a long list of parameters but most of them are self-explanatory
     * \c temp3 = units / one real pixel
     * \c temp4 = units / one imaginary pixel
     */
    // init r g and b
    uint r = 0, g = 0, b = 0;
    // initialize variables
    mpfr_t px, py, temp, temp2, re, im, wre, wim;
    mpfr_init2(px, lowerprec); mpfr_init2(py, lowerprec); 
    mpfr_init2(temp, precision); mpfr_init2(temp2, precision);
    mpfr_init2(re, precision); mpfr_init2(im, precision); 
    mpfr_init2(wre, precision); mpfr_init2(wim, precision);
    // this is the loop itself
    for (uint x = start_x; x < end_x; ++x) {
        mpfr_set_d(px, x, MPFR_RNDN); // a mpfr_t x value
        for (uint y = 0; y < height; ++y) {
            mpfr_set_d(py, y, MPFR_RNDN); // a mpfr_t y value
            mpfr_mul(temp, temp3, px, MPFR_RNDN); mpfr_add(re, minre, temp, MPFR_RNDN); // calculate re
            mpfr_mul(temp, temp4, py, MPFR_RNDN); mpfr_add(im, minim, temp, MPFR_RNDN); // calculate im
            mpfr_set_d(wre, 0, MPFR_RNDN); mpfr_set_d(wim, 0, MPFR_RNDN);
            uint t = 0; // the time
            while (t < max_iterations) {
                // Burning Ship fractal mod
                mpfr_abs(wim, wim, MPFR_RNDN);
                // do the real calc
                mpfr_sqr(temp, wre, MPFR_RNDN); mpfr_sqr(temp2, wim, MPFR_RNDN);
                mpfr_sub(temp, temp, temp2, MPFR_RNDN);
                mpfr_add(temp, temp, re, MPFR_RNDN); // add term
                // do the imaginary calc
                mpfr_mul(temp2, wim, wre, MPFR_RNDN); mpfr_mul(temp2, temp2, two, MPFR_RNDN);
                mpfr_add(temp2, temp2, im, MPFR_RNDN); // add term
                mpfr_set(wre, temp, MPFR_RNDN); mpfr_set(wim, temp2, MPFR_RNDN); // set values
                // calculate bailout
                mpfr_mul(temp, wre, wre, MPFR_RNDN); mpfr_mul(temp2, wim, wim, MPFR_RNDN);
                mpfr_add(temp, temp, temp2, MPFR_RNDN);
                if (mpfr_cmp(temp, bailout) > 0) break;
                ++t;
            }
            // store pixels
            if (t == max_iterations) {
                std::lock_guard<std::mutex> lock(image_mutex);
                image[height - y - 1][x] = png::rgb_pixel(0, 0, 0);
                continue;
            }
            hsvToRgb(((double) t) / colorlength, 1, 1, r, g, b);
            {
                std::lock_guard<std::mutex> lock(image_mutex);
                image[height - y - 1][x] = png::rgb_pixel(r, g, b);
            }
        }
    }
    // clear the (small amount of) mpfr_t variables we defined
    mpfr_clear(px); mpfr_clear(py); mpfr_clear(temp); mpfr_clear(temp2);
    mpfr_clear(re); mpfr_clear(im);
    mpfr_clear(wre); mpfr_clear(wim);
}

int mandelbrot(uint height, uint width,
    u_long precision, u_long loprec, u_long lowerprec,
    const char * CENTERRE, const char * CENTERIM,
    const char * BAILOUT, const char * SCALE,
    uint max_iterations, uint colorlength
){
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
    mpfr_t centerre, centerim, minre, dre, minim, dim;
    mpfr_t bailout, one, two, temp3, temp4, aspect, scale;
    mpfr_t lotemp, pwidth, pheight;
    char *endptr = NULL;
    // init
    mpfr_init2(centerre, precision); mpfr_init2(centerim, precision);
    mpfr_init2(aspect, loprec); mpfr_init2(scale, loprec);
    mpfr_init2(lotemp, loprec); mpfr_init2(bailout, loprec);
    mpfr_init2(one, lowerprec); mpfr_init2(two, lowerprec);
    mpfr_init2(temp3, precision); mpfr_init2(temp4, precision);
    mpfr_init2(minre, precision); mpfr_init2(dre, precision);
    mpfr_init2(minim, precision); mpfr_init2(dim, precision);
    mpfr_init2(pwidth, lowerprec); mpfr_init2(pheight, lowerprec);
    // values
    mpfr_set_d(pwidth, width, MPFR_RNDN); mpfr_set_d(pheight, height, MPFR_RNDN);
    // string values
    mpfr_strtofr(centerre, CENTERRE, &endptr, 10, MPFR_RNDN);
    mpfr_strtofr(centerim, CENTERIM, &endptr, 10, MPFR_RNDN);
    mpfr_strtofr(bailout, BAILOUT, &endptr, 10, MPFR_RNDN);
    mpfr_strtofr(scale, SCALE, &endptr, 10, MPFR_RNDN);
    // aspect ratio needs a calculation
    mpfr_set_d(aspect, width / height, MPFR_RNDN);
    // it's useful for us to have integer values here
    mpfr_set_d(one, 1, MPFR_RNDN); mpfr_set_d(two, 2, MPFR_RNDN);
    // calculation it looks like machine code :(
    // calculate real range
    mpfr_div(lotemp, aspect, scale, MPFR_RNDN);
    mpfr_sub(minre, centerre, lotemp, MPFR_RNDN);
    mpfr_mul(dre, two, lotemp, MPFR_RNDN);
    // calculation imaginary range
    mpfr_div(lotemp, one, scale, MPFR_RNDN);
    mpfr_sub(minim, centerim, lotemp, MPFR_RNDN);
    mpfr_mul(dim, two, lotemp, MPFR_RNDN);
    // calculate temp3 and temp4 for speed
    mpfr_div(temp3, dre, pwidth, MPFR_RNDN);
    mpfr_div(temp4, dim, pheight, MPFR_RNDN);
    png::image<png::rgb_pixel> image(width, height);
    std::cout
        << "CONFIG: ceneterre:" << CENTERRE
        << "centerim:" << CENTERIM
        << "scale:" << SCALE
        << "precision:" << precision
        << std::endl;
    // threads work
    uint num_threads = std::thread::hardware_concurrency();
    uint chunk_size = width / num_threads;
    std::vector<std::thread> threads;
    for (uint i = 0; i < num_threads; ++i) {
        uint start_x = i * chunk_size;
        uint end_x = (i == num_threads - 1) ? width : (i + 1) * chunk_size;
        threads.push_back(std::thread(mandelbrot_compute_pixels, start_x, end_x, height, width, 
                                      max_iterations, colorlength, precision, lowerprec,
                                      minre, dre, minim, dim,
                                      scale, std::ref(image),
                                      one, two, bailout, temp3, temp4));
    }
    // Join all threads
    std::for_each(threads.begin(), threads.end(), [](std::thread &t) { t.join(); });
    try {
        image.write("img.png");
        std::cout << "Image 'img.png' created successfully." << std::endl;
    } catch (const png::error& e) {
        std::cerr << "Error writing PNG: " << e.what() << std::endl;
        return 1;
    }
    // clear the (relatively more) variables we used
    mpfr_clear(centerre); mpfr_clear(centerim); mpfr_clear(aspect); mpfr_clear(scale);
    mpfr_clear(lotemp); mpfr_clear(bailout);
    mpfr_clear(one); mpfr_clear(two); mpfr_clear(temp3); mpfr_clear(temp4);
    mpfr_clear(minre); mpfr_clear(dre); mpfr_clear(minim); mpfr_clear(dim);
    mpfr_clear(pwidth); mpfr_clear(pheight);
    return 0;
}

int burningship(uint height, uint width,
    u_long precision, u_long loprec, u_long lowerprec,
    const char * CENTERRE, const char * CENTERIM,
    const char * BAILOUT, const char * SCALE,
    uint max_iterations, uint colorlength
){
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
    mpfr_t centerre, centerim, minre, dre, minim, dim;
    mpfr_t bailout, one, two, temp3, temp4, aspect, scale;
    mpfr_t lotemp, pwidth, pheight;
    char *endptr = NULL;
    // init
    mpfr_init2(centerre, precision); mpfr_init2(centerim, precision);
    mpfr_init2(aspect, loprec); mpfr_init2(scale, loprec);
    mpfr_init2(lotemp, loprec); mpfr_init2(bailout, loprec);
    mpfr_init2(one, lowerprec); mpfr_init2(two, lowerprec);
    mpfr_init2(temp3, precision); mpfr_init2(temp4, precision);
    mpfr_init2(minre, precision); mpfr_init2(dre, precision);
    mpfr_init2(minim, precision); mpfr_init2(dim, precision);
    mpfr_init2(pwidth, lowerprec); mpfr_init2(pheight, lowerprec);
    // values
    mpfr_set_d(pwidth, width, MPFR_RNDN); mpfr_set_d(pheight, height, MPFR_RNDN);
    // string values
    mpfr_strtofr(centerre, CENTERRE, &endptr, 10, MPFR_RNDN);
    mpfr_strtofr(centerim, CENTERIM, &endptr, 10, MPFR_RNDN);
    mpfr_strtofr(bailout, BAILOUT, &endptr, 10, MPFR_RNDN);
    mpfr_strtofr(scale, SCALE, &endptr, 10, MPFR_RNDN);
    // aspect ratio needs a calculation
    mpfr_set_d(aspect, width / height, MPFR_RNDN);
    // it's useful for us to have integer values here
    mpfr_set_d(one, 1, MPFR_RNDN); mpfr_set_d(two, 2, MPFR_RNDN);
    // calculation it looks like machine code :(
    // calculate real range
    mpfr_div(lotemp, aspect, scale, MPFR_RNDN);
    mpfr_sub(minre, centerre, lotemp, MPFR_RNDN);
    mpfr_mul(dre, two, lotemp, MPFR_RNDN);
    // calculation imaginary range
    mpfr_div(lotemp, one, scale, MPFR_RNDN);
    mpfr_sub(minim, centerim, lotemp, MPFR_RNDN);
    mpfr_mul(dim, two, lotemp, MPFR_RNDN);
    // calculate temp3 and temp4 for speed
    mpfr_div(temp3, dre, pwidth, MPFR_RNDN);
    mpfr_div(temp4, dim, pheight, MPFR_RNDN);
    png::image<png::rgb_pixel> image(width, height);
    std::cout
        << "CONFIG: ceneterre:" << CENTERRE
        << "centerim:" << CENTERIM
        << "scale:" << SCALE
        << "precision:" << precision
        << std::endl;
    // threads work
    uint num_threads = std::thread::hardware_concurrency();
    uint chunk_size = width / num_threads;
    std::vector<std::thread> threads;
    for (uint i = 0; i < num_threads; ++i) {
        uint start_x = i * chunk_size;
        uint end_x = (i == num_threads - 1) ? width : (i + 1) * chunk_size;
        threads.push_back(std::thread(mandelbrot_compute_pixels, start_x, end_x, height, width, 
                                      max_iterations, colorlength, precision, lowerprec,
                                      minre, dre, minim, dim,
                                      scale, std::ref(image),
                                      one, two, bailout, temp3, temp4));
    }
    // Join all threads
    std::for_each(threads.begin(), threads.end(), [](std::thread &t) { t.join(); });
    try {
        image.write("img.png");
        std::cout << "Image 'img.png' created successfully." << std::endl;
    } catch (const png::error& e) {
        std::cerr << "Error writing PNG: " << e.what() << std::endl;
        return 1;
    }
    // clear the (relatively more) variables we used
    mpfr_clear(centerre); mpfr_clear(centerim); mpfr_clear(aspect); mpfr_clear(scale);
    mpfr_clear(lotemp); mpfr_clear(bailout);
    mpfr_clear(one); mpfr_clear(two); mpfr_clear(temp3); mpfr_clear(temp4);
    mpfr_clear(minre); mpfr_clear(dre); mpfr_clear(minim); mpfr_clear(dim);
    mpfr_clear(pwidth); mpfr_clear(pheight);
    return 0;
}
