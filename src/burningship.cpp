// mainmandelbrot: the file for the burningship fractal
// Copyright (C) 2025-26  Edward Han
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
// along with this program.  If not, see <https://www.gnu.org/licenses/>
#define BURNINGSHIP
#include "common.h"
#include "util.h"
#ifdef BURNINGSHIP
double burningship_iters(
#else
double mandelbrot_iters(
#endif
    double re, double im, uint64_t max_iters, double bailout) {
    uint64_t t = 0;
    double temp, temp2, wre = 0, wim = 0, wre2 = 0, wim2 = 0, im2 = im * im;
#ifndef BURNINGSHIP
    temp = re - 0.25;
    double q = std::fma(temp, temp, im2);
    temp2 = q * (q + temp);
    temp2 = temp2 + temp2;
    if (temp2 + temp2 <= im2) {
        return max_iters;
    }
    temp += 1.25;
    if (std::fma(temp, temp, im2) <= 0.0625) {
        return max_iters;
    }
#endif
#ifdef TESTMASK
    return max_iters - 1;
#endif
    do {
        wre2 = wre * wre;
        wim2 = wim * wim;
        temp = wim * wre;
        wim = temp + temp + im;
        wre = wre2 - wim2 + re;
#ifdef BURNINGSHIP
        wim = std::abs(wim);
        wre = std::abs(wre);
#endif
        ++t;
    } while (t < max_iters && wim2 + wre2 < bailout);
    return t == max_iters ? t : (double)t - log2(log(wre2 + wim2));
}

#ifdef BURNINGSHIP
double burningship_iters(
#else
double mandelbrot_iters(
#endif
    mpfr_t re, mpfr_t im, uint64_t max_iters, mpfr_t bailout,
    u_long precision) {
    uint64_t t = 0;
    mpfr_t temp, temp2, wre, wim, wre2, wim2, im2;
    mpfr_inits2(precision, temp, temp2, wre, wim, wre2, wim2, im2, nullptr);
    mpfr_set_d(wre, 0, MPFR_RNDN);
    mpfr_set_d(wim, 0, MPFR_RNDN);
    mpfr_set_d(wre2, 0, MPFR_RNDN);
    mpfr_set_d(wim2, 0, MPFR_RNDN);
    mpfr_mul(im2, im, im, MPFR_RNDN);
#ifndef BURNINGSHIP
    mpfr_sub_d(temp, re, 0.25, MPFR_RNDN);
    mpfr_t q;
    mpfr_init2(q, precision);
    mpfr_fma(q, temp, temp, im2, MPFR_RNDN);
    mpfr_add(temp2, q, temp, MPFR_RNDN);
    mpfr_mul(temp2, q, temp2, MPFR_RNDN);
    mpfr_add(temp2, temp2, temp2, MPFR_RNDN);
    mpfr_add(temp2, temp2, temp2, MPFR_RNDN);
    if (mpfr_lessequal_p(temp2, im2)) {
        return max_iters;
    }
    mpfr_add_d(temp, temp, 1.25, MPFR_RNDN);
    mpfr_fma(temp2, temp, temp, im2, MPFR_RNDN);
    mpfr_set_d(temp, 0.0625, MPFR_RNDN);
    if (mpfr_lessequal_p(temp2, temp)) {
        return max_iters;
    }
#endif
    do {
        mpfr_sqr(wre2, wre, MPFR_RNDN);
        mpfr_sqr(wim2, wim, MPFR_RNDN);
        mpfr_mul(temp, wim, wre, MPFR_RNDN);
        mpfr_add(wim, temp, temp, MPFR_RNDN);
        mpfr_add(wim, wim, im, MPFR_RNDN);
        mpfr_sub(wre, wre2, wim2, MPFR_RNDN);
        mpfr_add(wre, wre, re, MPFR_RNDN);
#ifdef BURNINGSHIP
        mpfr_abs(wre, wre, MPFR_RNDN);
        mpfr_abs(wim, wim, MPFR_RNDN);
#endif
        mpfr_add(temp, wre2, wim2, MPFR_RNDN);
        ++t;
    } while (t < max_iters && mpfr_less_p(temp, bailout));
    double result =
        t == max_iters ? t : (double)t - log2(log(mpfr_get_d(temp, MPFR_RNDN)));
    mpfr_clears(temp, temp2, wre, wim, wre2, wim2, im2, nullptr);
    return result;
}

#ifdef BURNINGSHIP
void burningship_work_double(std::atomic<uint64_t> &x, uint64_t height,
#else
void mandelbrot_work_double(std::atomic<uint64_t> &x, uint64_t height,
#endif
                             uint64_t width, uint64_t max_iterations,
                             uint64_t colorlength, double minre, double minim,
                             png::image<png::rgb_pixel> &image, double bailout,
                             double temp3, double temp4) {
    uint64_t size = width * height;
    uint64_t work_x = x.fetch_add(1);
    while (work_x < size) {
        int64_t x_ = work_x % width, y_ = work_x / width;
        double re = std::fma(temp3, x_, minre), im = std::fma(temp4, y_, minim);
        double iters =
#ifdef BURNINGSHIP
            burningship_iters(
#else
            mandelbrot_iters(
#endif
                re, im, max_iterations, bailout);
        // store pixels
        if (iters == (double)max_iterations) {
            image[height - y_ - 1][x_] = png::rgb_pixel(0, 0, 0);
        } else {
            double tp = std::fmod(iters / colorlength, 1.);
            image[height - y_ - 1][x_] = grad_0(tp);
        }
        work_x = x.fetch_add(1);
    }
}

#ifdef BURNINGSHIP
void burningship_work(
#else
void mandelbrot_work(
#endif
    std::atomic<uint64_t> &x, uint64_t height, uint64_t width,
    uint64_t max_iterations, uint64_t colorlength, mpfr_t minre, mpfr_t dre,
    mpfr_t minim, mpfr_t dim, mpfr_t scale, png::image<png::rgb_pixel> &image,
    mpfr_t bailout, mpfr_t temp3, mpfr_t temp4, u_long precision) {
    uint64_t size = width * height;
    uint64_t work_x = x.fetch_add(1);
    mpfr_t re, im;
    mpfr_inits2(precision, re, im, nullptr);
    while (work_x < size) {
        int64_t x_ = work_x % width, y = work_x / width;
        mpfr_set_d(re, x_, MPFR_RNDN);
        mpfr_fma(re, temp3, re, minre, MPFR_RNDN);
        mpfr_set_d(im, y, MPFR_RNDN);
        mpfr_fma(im, temp4, im, minim, MPFR_RNDN);
        double iters =
#ifdef BURNINGSHIP
            burningship_iters(
#else
            mandelbrot_iters(
#endif
                re, im, max_iterations, bailout, precision);
        if (iters == (double)max_iterations) {
            image[height - y - 1][x_] = png::rgb_pixel(0, 0, 0);
        } else {
            double tp = std::fmod(iters / colorlength, 1.);
            image[height - y - 1][x_] = grad_0(tp);
        }
        work_x = x.fetch_add(1);
    }
    mpfr_clears(re, im, nullptr);
}
#ifdef BURNINGSHIP
int burningship(uint64_t height, uint64_t width, u_long precision,
                u_long loprec,
#else
int mandelbrot(uint64_t height, uint64_t width, u_long precision,
#endif
                const char *CENTERRE, const char *CENTERIM, const char *BAILOUT,
                const char *SCALE, uint64_t max_iterations,
                uint64_t colorlength, std::string filename) {
    mpfr_t centerre, centerim, minre, dre, minim, dim;
    mpfr_t bailout, temp3, temp4, aspect, scale;
    mpfr_t lotemp, pwidth, pheight;
    // init
    mpfr_inits2(precision, centerre, centerim, temp3, temp4, minre, dre, minim,
                dim, nullptr);
    mpfr_inits2(27, aspect, scale, lotemp, bailout, pwidth, pheight, nullptr);
    // values
    mpfr_set_d(pwidth, width, MPFR_RNDN);
    mpfr_set_d(pheight, height, MPFR_RNDN);
    // string values
    mpfr_strtofr(centerre, CENTERRE, nullptr, 10, MPFR_RNDN);
    mpfr_strtofr(centerim, CENTERIM, nullptr, 10, MPFR_RNDN);
    mpfr_strtofr(bailout, BAILOUT, nullptr, 10, MPFR_RNDN);
    mpfr_strtofr(scale, SCALE, nullptr, 10, MPFR_RNDN);
    // aspect ratio needs a calculation
    mpfr_set_d(aspect, width / (double)height, MPFR_RNDN);
    // it's useful for us to have integer values here
    // calculation it looks like machine code :(
    // calculate real range
    mpfr_div(lotemp, aspect, scale, MPFR_RNDN);
    mpfr_sub(minre, centerre, lotemp, MPFR_RNDN);
    mpfr_add(dre, lotemp, lotemp, MPFR_RNDN);
    // calculation imaginary range
    mpfr_ui_div(lotemp, 1, scale, MPFR_RNDN);
    mpfr_sub(minim, centerim, lotemp, MPFR_RNDN);
    mpfr_add(dim, lotemp, lotemp, MPFR_RNDN);
    // calculate temp3 and temp4 for speed
    mpfr_div(temp3, dre, pwidth, MPFR_RNDN);
    mpfr_div(temp4, dim, pheight, MPFR_RNDN);
    png::image<png::rgb_pixel> image(width, height);
    // std::cout << "CONFIG: ceneterre:" << CENTERRE << "centerim:" << CENTERIM
    //           << "scale:" << SCALE << "precision:" << precision << std::endl;
    // threads work
    uint64_t num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    std::atomic<uint64_t> x{0};
    for (uint64_t i = 0; i < num_threads; ++i) {
        threads.push_back(std::thread(
#ifdef BURNINGSHIP
            burningship_work, std::ref(x), height, width, max_iterations,
#else
            mandelbrot_work, std::ref(x), height, width, max_iterations,
#endif
            colorlength, minre, dre, minim, dim, scale, std::ref(image),
            bailout, temp3, temp4, precision));
    }
    // Join all threads
    std::for_each(threads.begin(), threads.end(),
                  [](std::thread &t) { t.join(); });
    try {
        image.write(filename);
    } catch (const png::error &e) {
        std::cerr << "Error writing PNG: " << e.what() << std::endl;
        return 1;
    }
    // clear the (relatively more) variables we used
    mpfr_clears(centerre, centerim, temp3, temp4, minre, dre, minim, dim,
                aspect, scale, lotemp, bailout, pwidth, pheight, nullptr);
    return 0;
}

#ifdef BURNINGSHIP
int burningship(uint64_t height, uint64_t width, double centerre,
#else
int mandelbrot(uint64_t height, uint64_t width, double centerre,
#endif
                double centerim, double bailout, double scale,
                uint64_t max_iterations, uint64_t colorlength,
                std::string filename) {
    double minre, dre, minim, dim, temp3, temp4, aspect, lotemp;
    aspect = width / (double)height;
    // calculate real
    lotemp = aspect / scale;
    minre = centerre - lotemp;
    dre = lotemp + lotemp;
    // calculation imaginary range
    lotemp = 1 / scale;
    minim = centerim - lotemp;
    dim = lotemp + lotemp;
    // calculate temp3 and temp4 for speed (per pixel values)
    temp3 = dre / width;
    temp4 = dim / height;
    png::image<png::rgb_pixel> image(width, height);
    // std::cout << "CONFIG: ceneterre:" << centerre << "centerim:" << centerim
    //   << "scale:" << scale << "(double)" << std::endl;
    // threads work
    uint64_t num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    std::atomic<uint64_t> x{0};

    for (uint64_t i = 0; i < num_threads; ++i) {
        threads.push_back(
#ifdef BURNINGSHIP
            std::thread(burningship_work_double, std::ref(x), height, width,
#else
            std::thread(mandelbrot_work_double, std::ref(x), height, width,
#endif
                        max_iterations, colorlength, minre, minim,
                        std::ref(image), bailout, temp3, temp4));
    }
    // Join all threads
    std::for_each(threads.begin(), threads.end(),
                  [](std::thread &t) { t.join(); });
    try {
        image.write(filename);
    } catch (const png::error &e) {
        std::cerr << "Error writing PNG: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
