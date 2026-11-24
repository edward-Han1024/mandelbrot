#include "common.h"
#pragma once

class Pchip {
 public:
  Pchip(std::vector<double> x, std::vector<double> y);
  Pchip();
  double get_value(double x);
  double operator()(double x);
 private:
  static int sign(double value);
  static double endpoint_slope(double h0, double h1, double delta0,
                               double delta1);
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> slopes_;
};
