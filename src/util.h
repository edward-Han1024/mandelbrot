#include "common.h"
#pragma once

/**
 * @brief Convertes hsv value into an rgb value.
 *
 * @param h (double) The hue value. The range is [0, 1), but it can loop
 * around.
 * @param s (double) The saturation of the color. The range is [0, 1)
 * @param v (double) The value of the color. The range is [0, 1)
 * @param r (uint&) A pointer to the uint to store the value of r. The range
 * will be [0, 255].
 * @param g (uint&) A pointer to the uint to store the value of g. The range
 * will be [0, 255].
 * @param b (uint&) A pointer to the uint to store the value of b. The range
 * will be [0, 255].
 */
void hsvToRgb(double h, double s, double v, uint &r, uint &g, uint &b);
png::rgb_pixel grad_0(double t);
// Grayscale.
png::rgb_pixel grad_1(double t);
// Viridis.
png::rgb_pixel grad_2(double t);
// Inferno.
png::rgb_pixel grad_3(double t);
// Coolwarm.
png::rgb_pixel grad_4(double t);
png::rgb_pixel grad_5(double t);
png::rgb_pixel grad_6(double t);
png::rgb_pixel grad_7(double t);
png::rgb_pixel grad_8(double t);
png::rgb_pixel grad_9(double t);