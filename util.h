#include "common.h"
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
void
hsvToRgb (double h, double s, double v, uint &r, uint &g, uint &b);
