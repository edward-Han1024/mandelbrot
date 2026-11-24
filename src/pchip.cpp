#include "pchip.h"
#include "common.h"

Pchip::Pchip(std::vector<double> x, std::vector<double> y)
    : x_(std::move(x)), y_(std::move(y)), slopes_(x_.size()) {
    const std::size_t n = x_.size();
    std::vector<double> h(n - 1);
    std::vector<double> secant(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x_[i + 1] - x_[i];
        secant[i] = (y_[i + 1] - y_[i]) / h[i];
    }

    if (n == 2) {
        slopes_[0] = secant[0];
        slopes_[1] = secant[0];
        return;
    }

    slopes_[0] = endpoint_slope(h[0], h[1], secant[0], secant[1]);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        const double previous = secant[i - 1];
        const double next = secant[i];
        if (previous == 0.0 || next == 0.0 ||
            ((previous > 0.0) != (next > 0.0))) {
            slopes_[i] = 0.0;
        } else {
            const double w1 = 2.0 * h[i] + h[i - 1];
            const double w2 = h[i] + 2.0 * h[i - 1];
            slopes_[i] = (w1 + w2) / (w1 / previous + w2 / next);
        }
    }
    slopes_[n - 1] =
        endpoint_slope(h[n - 2], h[n - 3], secant[n - 2], secant[n - 3]);
}

double Pchip::get_value(double x) {
    const std::size_t n = x_.size();
    std::size_t i;
    if (x <= x_.front()) {
        i = 0;
    } else if (x >= x_.back()) {
        i = n - 2;
    } else {
        i = static_cast<std::size_t>(std::upper_bound(x_.begin(), x_.end(), x) -
                                     x_.begin() - 1);
    }
    
    const double h = x_[i + 1] - x_[i];
    const double t = (x - x_[i]) / h;
    const double t2 = t * t;
    const double t3 = t2 * t;
    return (2.0 * t3 - 3.0 * t2 + 1.0) * y_[i] +
           (t3 - 2.0 * t2 + t) * h * slopes_[i] +
           (-2.0 * t3 + 3.0 * t2) * y_[i + 1] + (t3 - t2) * h * slopes_[i + 1];
}

int Pchip::sign(double value) { return (value > 0.0) - (value < 0.0); }

double Pchip::endpoint_slope(double h0, double h1, double delta0,
                             double delta1) {
    double slope = ((2.0 * h0 + h1) * delta0 - h0 * delta1) / (h0 + h1);
    if (sign(slope) != sign(delta0)) {
        return 0.0;
    }
    if (sign(delta0) != sign(delta1) &&
        std::abs(slope) > std::abs(3.0 * delta0)) {
        return 3.0 * delta0;
    }
    return slope;
}

double Pchip::operator()(double x) { return get_value(x); }

Pchip::Pchip() : Pchip({0, 0.5, 1}, {0, 0.5, 1}) {}

