// mainmandelbrot: the file describing the mandelbrot fractal
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

#include "common.h"
#include "util.h"

std::mutex mandel_image_mutex;
std::mutex mandel_work_mutex;

void
mandelbrot_compute_pixel (uint x, uint y, uint height, uint width,
                          uint max_iterations, uint colorlength,
                          uint precision, uint lowerprec, mpfr_t minre,
                          mpfr_t dre, mpfr_t minim, mpfr_t dim, mpfr_t scale,
                          png::image<png::rgb_pixel> &image, mpfr_t bailout,
                          mpfr_t temp3, mpfr_t temp4)
{
  if (y == height)
    {
      return;
    }
  // std::cout << "Called on x: " << x << "and y: " << y << std::endl;
  // init r g and b
  uint64_t r = 0, g = 0, b = 0;
  // initialize variables
  mpfr_t px, py, temp, temp2, re, im, wre, wim;
  mpfr_inits2 (lowerprec, px, py, nullptr);
  mpfr_inits2 (precision, temp, temp2, re, im, wre, wim, nullptr);
  // this is the loop itself
  mpfr_set_d (px, x, MPFR_RNDN); // a mpfr_t x value
  mpfr_set_d (py, y, MPFR_RNDN); // a mpfr_t y value
  mpfr_mul (temp, temp3, px, MPFR_RNDN);
  mpfr_add (re, minre, temp, MPFR_RNDN); // calculate re
  mpfr_mul (temp, temp4, py, MPFR_RNDN);
  mpfr_add (im, minim, temp, MPFR_RNDN); // calculate im
  mpfr_set_d (wre, 0, MPFR_RNDN);
  mpfr_set_d (wim, 0, MPFR_RNDN);
  uint64_t t = 0; // the time
  while (t < max_iterations)
    {
      // do the real calc
      mpfr_sqr (temp, wre, MPFR_RNDN);
      mpfr_sqr (temp2, wim, MPFR_RNDN);
      mpfr_sub (temp, temp, temp2, MPFR_RNDN);
      mpfr_add (temp, temp, re, MPFR_RNDN); // add term
      // do the imaginary calc
      mpfr_mul (temp2, wim, wre, MPFR_RNDN);
      mpfr_add (temp2, temp2, temp2, MPFR_RNDN);
      mpfr_add (temp2, temp2, im, MPFR_RNDN); // add term
      mpfr_set (wre, temp, MPFR_RNDN);
      mpfr_set (wim, temp2, MPFR_RNDN); // set values
      // calculate bailout
      mpfr_sqr (temp, wre, MPFR_RNDN);
      mpfr_sqr (temp2, wim, MPFR_RNDN);
      mpfr_add (temp, temp, temp2, MPFR_RNDN);
      if (mpfr_cmp (temp, bailout) > 0)
        break;
      ++t;
    }
  // store pixels
  if (t == max_iterations)
    {
      std::lock_guard<std::mutex> lock (mandel_image_mutex);
      image[height - y - 1][x] = png::rgb_pixel (0, 0, 0);
      return;
    }
  hsvToRgb (((double)t) / colorlength, 1, 1, r, g, b);
  {
    std::lock_guard<std::mutex> lock (mandel_image_mutex);
    image[height - y - 1][x] = png::rgb_pixel (r, g, b);
  }
  // clear the (small amount of) mpfr_t variables we defined
  mpfr_clears (px, py, temp, temp2, re, im, wre, wim, nullptr);
}

void
mandelbrot_compute_pixel (uint64_t x, uint64_t y, uint64_t height,
                          uint64_t width, uint64_t max_iterations,
                          uint64_t colorlength, double minre, double dre,
                          double minim, double dim, double scale,
                          png::image<png::rgb_pixel> &image, double bailout,
                          double temp3, double temp4)
{
  if (y == height)
    return;
  //  std::cout << "Called on x: " << x << "and y: " << y << std::endl;
  // init r g and b
  uint64_t r = 0, g = 0, b = 0;
  // initialize variables
  double temp, temp2, re, im, wre, wim;
  // this is the loop itself
  re = minre + temp3 * x; // calculate re
  im = minim + temp4 * y; // calculate im
  wre = 0;
  wim = 0;
  uint64_t t = 0; // the time
  // do bulb checking.
  double q = (re - 0.25) * (re - 0.25) + y * y;
  if (q * (q + re - 0.25) <= 0.25 * y * y)
    {
      std::lock_guard<std::mutex> lock (mandel_image_mutex);
      image[height - y - 1][x] = png::rgb_pixel (0, 0, 0);
      return;
    } // cardiod
  if ((re + 1) * (re + 1) + y * y <= 0.0625)
    {
      std::lock_guard<std::mutex> lock (mandel_image_mutex);
      image[height - y - 1][x] = png::rgb_pixel (0, 0, 0);
      return;
    } // circle
  while (t < max_iterations)
    {
      // do the real calc
      temp = wre * wre - wim * wim + re;
      // do the imaginary calc
      temp2 = wim * wre;
      temp2 = temp2 + temp2 + im;
      wre = temp;
      wim = temp2;
      // calculate bailout
      temp = wre * wre + wim * wim;
      if (temp > bailout)
        break;
      ++t;
    }
  // store pixels
  if (t == max_iterations)
    {
      std::lock_guard<std::mutex> lock (mandel_image_mutex);
      image[height - y - 1][x] = png::rgb_pixel (0, 0, 0);
      return;
    }
  hsvToRgb (((double)t) / colorlength, 1, 1, r, g, b);
  {
    std::lock_guard<std::mutex> lock (mandel_image_mutex);
    image[height - y - 1][x] = png::rgb_pixel (r, g, b);
  }
}

void
mandelbrot_work (uint64_t &x, uint64_t &y, uint64_t height, uint64_t width,
                 uint64_t max_iterations, uint64_t colorlength,
                 uint64_t precision, uint64_t lowerprec, mpfr_t minre,
                 mpfr_t dre, mpfr_t minim, mpfr_t dim, mpfr_t scale,
                 png::image<png::rgb_pixel> &image, mpfr_t bailout,
                 mpfr_t temp3, mpfr_t temp4)
{
  while ((x < width) && (y < height))
    {
      // std::cout << "[" << x << y << "]" <<  std::endl;
      uint64_t work_x, work_y;
      { // Scope for work update.
        std::lock_guard<std::mutex> lock (mandel_work_mutex);
        work_x = x;
        work_y = y;
        x++;
        if (x >= width)
          {
            x = 0;
            y++;
          }
      }
      mandelbrot_compute_pixel (work_x, work_y, height, width, max_iterations,
                                colorlength, precision, lowerprec, minre, dre,
                                minim, dim, scale, image, bailout, temp3,
                                temp4);
    }
}

void
mandelbrot_work_double (uint64_t &x, uint64_t &y, uint64_t height,
                        uint64_t width, uint64_t max_iterations,
                        uint64_t colorlength, double minre, double dre,
                        double minim, double dim, double scale,
                        png::image<png::rgb_pixel> &image, double bailout,
                        double temp3, double temp4)
{
  while ((x < width) && (y < height))
    {
      // std::cout << "[" << x << y << "]" <<  std::endl;
      uint64_t work_x, work_y;
      { // Scope for work update.
        std::lock_guard<std::mutex> lock (mandel_work_mutex);
        work_x = x;
        work_y = y;
        x++;
        if (x == width)
          {
            x = 0;
            y++;
          }
      }
      mandelbrot_compute_pixel (work_x, work_y, height, width, max_iterations,
                                colorlength, minre, dre, minim, dim, scale,
                                image, bailout, temp3, temp4);
    }
}

int
mandelbrot (uint64_t height, uint64_t width, u_long precision, u_long loprec,
            u_long lowerprec, const char *CENTERRE, const char *CENTERIM,
            const char *BAILOUT, const char *SCALE, uint64_t max_iterations,
            uint64_t colorlength)
{
  mpfr_t centerre, centerim, minre, dre, minim, dim;
  mpfr_t bailout, one, two, temp3, temp4, aspect, scale;
  mpfr_t lotemp, pwidth, pheight;
  // init
  mpfr_inits2 (precision, centerre, centerim, temp3, temp4, minre, dre, minim,
               dim, nullptr);
  mpfr_inits2 (loprec, aspect, scale, lotemp, bailout, pwidth, pheight,
               nullptr);
  mpfr_inits2 (lowerprec, one, two, nullptr);
  // values
  mpfr_set_d (pwidth, width, MPFR_RNDN);
  mpfr_set_d (pheight, height, MPFR_RNDN);
  // string values
  mpfr_strtofr (centerre, CENTERRE, nullptr, 10, MPFR_RNDN);
  mpfr_strtofr (centerim, CENTERIM, nullptr, 10, MPFR_RNDN);
  mpfr_strtofr (bailout, BAILOUT, nullptr, 10, MPFR_RNDN);
  mpfr_strtofr (scale, SCALE, nullptr, 10, MPFR_RNDN);
  // aspect ratio needs a calculation
  mpfr_set_d (aspect, width / (double)height, MPFR_RNDN);
  // it's useful for us to have integer values here
  // calculation it looks like machine code :(
  // calculate real range
  mpfr_div (lotemp, aspect, scale, MPFR_RNDN);
  mpfr_sub (minre, centerre, lotemp, MPFR_RNDN);
  mpfr_add (dre, lotemp, lotemp, MPFR_RNDN);
  // calculation imaginary range
  mpfr_ui_div (lotemp, 1, scale, MPFR_RNDN);
  mpfr_sub (minim, centerim, lotemp, MPFR_RNDN);
  mpfr_add (dim, lotemp, lotemp, MPFR_RNDN);
  // calculate temp3 and temp4 for speed
  mpfr_div (temp3, dre, pwidth, MPFR_RNDN);
  mpfr_div (temp4, dim, pheight, MPFR_RNDN);
  png::image<png::rgb_pixel> image (width, height);
  // std::cout << "CONFIG: ceneterre:" << CENTERRE << "centerim:" << CENTERIM
  //           << "scale:" << SCALE << "precision:" << precision << std::endl;
  // threads work
  uint64_t num_threads = std::thread::hardware_concurrency ();
  std::vector<std::thread> threads;
  uint64_t x = 0;
  uint64_t y = 0;

  for (uint64_t i = 0; i < num_threads; ++i)
    {
      threads.push_back (std::thread (
          mandelbrot_work, std::ref (x), std::ref (y), height, width,
          max_iterations, colorlength, precision, lowerprec, minre, dre, minim,
          dim, scale, std::ref (image), bailout, temp3, temp4));
    }
  // Join all threads
  std::for_each (threads.begin (), threads.end (),
                 [] (std::thread &t) { t.join (); });
  try
    {
      image.write ("img.png");
    }
  catch (const png::error &e)
    {
      std::cerr << "Error writing PNG: " << e.what () << std::endl;
      return 1;
    }
  // clear the (relatively more) variables we used
  mpfr_clears (centerre, centerim, aspect, scale, lotemp, bailout, one, two,
               temp4, temp4, minre, dre, minim, dim, pwidth, pheight, nullptr);
  return 0;
}

int
mandelbrot (uint64_t height, uint64_t width, double centerre, double centerim,
            double bailout, double scale, uint64_t max_iterations,
            uint64_t colorlength)
{
  double minre, dre, minim, dim, temp3, temp4, aspect, lotemp;
  // init
  // values
  aspect = width / (double)height;
  // calculation it looks like machine code :(
  // calculate real
  lotemp = aspect / scale;
  minre = centerre - lotemp;
  dre = lotemp + lotemp;
  // calculation imaginary range
  lotemp = 1 / scale;
  minim = centerim - lotemp;
  dim = 2 * lotemp;
  // calculate temp3 and temp4 for speed
  temp3 = dre / width;
  temp4 = dim / height;
  png::image<png::rgb_pixel> image (width, height);
  std::cout << "CONFIG: ceneterre:" << centerre << "centerim:" << centerim
            << "scale:" << scale << "(double)" << std::endl;
  // threads work
  uint64_t num_threads = std::thread::hardware_concurrency ();
  std::vector<std::thread> threads;
  uint64_t x = 0;
  uint64_t y = 0;

  for (uint64_t i = 0; i < num_threads; ++i)
    {
      threads.push_back (std::thread (
          mandelbrot_work_double, std::ref (x), std::ref (y), height, width,
          max_iterations, colorlength, minre, dre, minim, dim, scale,
          std::ref (image), bailout, temp3, temp4));
    }
  // Join all threads
  std::for_each (threads.begin (), threads.end (),
                 [] (std::thread &t) { t.join (); });
  try
    {
      image.write ("img.png");
    }
  catch (const png::error &e)
    {
      std::cerr << "Error writing PNG: " << e.what () << std::endl;
      return 1;
    }
  return 0;
}
