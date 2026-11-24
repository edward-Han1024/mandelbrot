#include "common.h"
#include "pchip.h"
#pragma once

class Gradient {
  public:
    Gradient(std::vector<double> c, std::vector<double> r,
             std::vector<double> g, std::vector<double> b);
    Gradient(std::string filename);
    png::rgb_pixel operator()(double t);
 private:
    Pchip r, g, b;
};
