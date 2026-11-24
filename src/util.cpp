#include "common.h"
#include "pchip.h"

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
    t -= std::floor(t);
    static std::vector<double> grad0_c{0.0, 0.16, 0.42, 0.6425, 0.8575, 1};
    static std::vector<double> grad0_r{0., 32., 237., 255., 0., 0.};
    static std::vector<double> grad0_g{7., 107., 255., 170., 2., 7.};
    static std::vector<double> grad0_b{100., 203., 255., 0., 0., 100.};
    static auto spline0_r = Pchip(grad0_c, grad0_r);
    static auto spline0_g = Pchip(grad0_c, grad0_g);
    static auto spline0_b = Pchip(grad0_c, grad0_b);
    return png::rgb_pixel(spline0_r.get_value(t), spline0_g.get_value(t),
                          spline0_b.get_value(t));
}

png::rgb_pixel grad_1(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad1_c{0.0, 0.5, 1.0};
    static const std::vector<double> grad1_rgb{0., 255., 0.};
    static Pchip spline1(grad1_c, grad1_rgb);
    const auto value = spline1.get_value(t);
    return png::rgb_pixel(value, value, value);
}

png::rgb_pixel grad_2(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad2_c{0.0, 0.25, 0.5, 0.75, 1.0};
    static const std::vector<double> grad2_r{68., 59., 33., 94., 68.};
    static const std::vector<double> grad2_g{1., 82., 145., 201., 1.};
    static const std::vector<double> grad2_b{84., 139., 140., 98., 84.};
    static Pchip spline2_r(grad2_c, grad2_r);
    static Pchip spline2_g(grad2_c, grad2_g);
    static Pchip spline2_b(grad2_c, grad2_b);
    return png::rgb_pixel(spline2_r.get_value(t), spline2_g.get_value(t),
                          spline2_b.get_value(t));
}

png::rgb_pixel grad_3(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad3_c{0.0, 0.25, 0.5, 0.75, 1.0};
    static const std::vector<double> grad3_r{0., 87., 188., 249., 0.};
    static const std::vector<double> grad3_g{0., 15., 55., 142., 0.};
    static const std::vector<double> grad3_b{4., 109., 84., 8., 4.};
    static Pchip spline3_r(grad3_c, grad3_r);
    static Pchip spline3_g(grad3_c, grad3_g);
    static Pchip spline3_b(grad3_c, grad3_b);
    return png::rgb_pixel(spline3_r.get_value(t), spline3_g.get_value(t),
                          spline3_b.get_value(t));
}

png::rgb_pixel grad_4(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad4_c{0.0, 0.25, 0.5, 0.75, 1.0};
    static const std::vector<double> grad4_r{59., 141., 221., 245., 59.};
    static const std::vector<double> grad4_g{76., 176., 221., 156., 76.};
    static const std::vector<double> grad4_b{192., 254., 221., 125., 192.};
    static Pchip spline4_r(grad4_c, grad4_r);
    static Pchip spline4_g(grad4_c, grad4_g);
    static Pchip spline4_b(grad4_c, grad4_b);
    return png::rgb_pixel(spline4_r.get_value(t), spline4_g.get_value(t),
                          spline4_b.get_value(t));
}

png::rgb_pixel grad_5(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad6_c{0.0, 0.3, 0.65, 0.85, 1.0};
    static const std::vector<double> grad6_r{0.0, 120.0, 255.0, 255.0, 255.0};
    static const std::vector<double> grad6_g{0.0, 10.0, 100.0, 230.0, 255.0};
    static const std::vector<double> grad6_b{0.0, 10.0, 0.0, 100.0, 255.0};
    static Pchip spline6_r(grad6_c, grad6_r);
    static Pchip spline6_g(grad6_c, grad6_g);
    static Pchip spline6_b(grad6_c, grad6_b);
    return png::rgb_pixel(spline6_r.get_value(t), spline6_g.get_value(t),
                          spline6_b.get_value(t));
}

png::rgb_pixel grad_6(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad10_c{0.0, 0.3, 0.6, 0.85, 1.0};
    static const std::vector<double> grad10_r{0.0, 0.0,   50.0,  200.0, 0.0};
    static const std::vector<double> grad10_g{10.0, 120.0, 220.0, 255.0, 10.0};
    static const std::vector<double> grad10_b{5.0,  50.0,  100.0, 200.0, 5.0};
    static Pchip spline10_r(grad10_c, grad10_r);
    static Pchip spline10_g(grad10_c, grad10_g);
    static Pchip spline10_b(grad10_c, grad10_b);
    return png::rgb_pixel(spline10_r.get_value(t), spline10_g.get_value(t),
                          spline10_b.get_value(t));
}

png::rgb_pixel grad_7(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad7_c{0.0, 0.166, 0.333, 0.5, 0.666, 0.833, 1.0};
    static const std::vector<double> grad7_r{255.0, 255.0, 0.0,   0.0,   0.0,   130.0, 255.0};
    static const std::vector<double> grad7_g{0.0,   127.0, 255.0, 255.0, 0.0,   0.0,   0.0};
    static const std::vector<double> grad7_b{0.0,   0.0,   0.0,   255.0, 255.0, 255.0, 0.0};
    static Pchip spline7_r(grad7_c, grad7_r);
    static Pchip spline7_g(grad7_c, grad7_g);
    static Pchip spline7_b(grad7_c, grad7_b);
    return png::rgb_pixel(spline7_r.get_value(t), spline7_g.get_value(t),
                          spline7_b.get_value(t));
}

png::rgb_pixel grad_8(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad8_c{0.0, 0.25, 0.5, 0.75, 1.0};
    static const std::vector<double> grad8_r{15.0, 130.0, 255.0, 0.0,   15.0};
    static const std::vector<double> grad8_g{5.0,  0.0,   20.0,  230.0, 5.0};
    static const std::vector<double> grad8_b{30.0, 180.0, 130.0, 255.0, 30.0};
    static Pchip spline8_r(grad8_c, grad8_r);
    static Pchip spline8_g(grad8_c, grad8_g);
    static Pchip spline8_b(grad8_c, grad8_b);
    return png::rgb_pixel(spline8_r.get_value(t), spline8_g.get_value(t),
                          spline8_b.get_value(t));
}

png::rgb_pixel grad_9(double t) {
    t -= std::floor(t);
    static const std::vector<double> grad9_c{0.0, 0.25, 0.5, 0.75, 1.0};
    static const std::vector<double> grad9_r{15.0, 120.0, 230.0, 250.0, 15.0};
    static const std::vector<double> grad9_g{10.0, 60.0,  125.0, 210.0, 10.0};
    static const std::vector<double> grad9_b{5.0,  20.0,  50.0,  140.0, 5.0};
    static Pchip spline9_r(grad9_c, grad9_r);
    static Pchip spline9_g(grad9_c, grad9_g);
    static Pchip spline9_b(grad9_c, grad9_b);
    return png::rgb_pixel(spline9_r.get_value(t), spline9_g.get_value(t),
                          spline9_b.get_value(t));
}


