#include "common.h"

void hsvToRgb(double h, double s, double v, uint &r, uint &g, uint &b) {
    double c = v * s;
    double hPrime = fmod(6.0 * h, 6.0); // Account for looping over
    double x = c * (1 - std::abs(fmod(hPrime, 2) - 1));
    double r1, g1, b1;
    // There are 6 cases for the conversion
    if (hPrime >= 0 && hPrime < 1) {
        r1 = c;
        g1 = x;
        b1 = 0;
    } else if (hPrime >= 1 && hPrime < 2) {
        r1 = x;
        g1 = c;
        b1 = 0;
    } else if (hPrime >= 2 && hPrime < 3) {
        r1 = 0;
        g1 = c;
        b1 = x;
    } else if (hPrime >= 3 && hPrime < 4) {
        r1 = 0;
        g1 = x;
        b1 = c;
    } else if (hPrime >= 4 && hPrime < 5) {
        r1 = x;
        g1 = 0;
        b1 = c;
    } else { // hPrime >= 5 && hPrime < 6
        r1 = c;
        g1 = 0;
        b1 = x;
    }
    // brightness correction
    double m = v - c;
    r = (r1 + m) * 255.0;
    g = (g1 + m) * 255.0;
    b = (b1 + m) * 255.0;
}

png::rgb_pixel grad_0(double t) {
    static std::vector<double> grad0_c{0.0, 0.16, 0.42, 0.6425, 0.8575, 1};
    static std::vector<double> grad0_r{0., 32., 237., 255., 0., 0.};
    static std::vector<double> grad0_g{7., 107., 255., 170., 2., 7.};
    static std::vector<double> grad0_b{100., 203., 255., 0., 0., 100.};
    static auto spline0_r = pchip(std::move(grad0_c), std::move(grad0_r));
    static std::vector<double> grad0_c1 = std::vector<double>{0.0, 0.16, 0.42, 0.6425, 0.8575, 1};
    static auto spline0_g = pchip(std::move(grad0_c1), std::move(grad0_g));
    static std::vector<double> grad0_c2 = std::vector<double>{0.0, 0.16, 0.42, 0.6425, 0.8575, 1};
    static auto spline0_b = pchip(std::move(grad0_c2), std::move(grad0_b));
    return png::rgb_pixel(spline0_r(t), spline0_g(t), spline0_b(t));
}