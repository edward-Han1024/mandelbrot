#include "common.h"

struct transform
{
  std::function<std::pair<double, double> (
      std::vector<std::pair<double, double>>, std::pair<double, double>)>
      transformation;
  std::vector<std::pair<double, double>> parameters;
  double weight;
  std::pair<double, double>
  apply_transform (std::pair<double, double> point)
  {
    return transformation (parameters, point);
  }
};

std::pair<double, double>
affine (std::vector<std::pair<double, double>> parameters,
        std::pair<double, double> point)
{
  /**
   * @brief the affine transformation aka the linear transformation
   */
  double X1 = parameters[0].first;
  double X2 = parameters[0].second;
  double Y1 = parameters[1].first;
  double Y2 = parameters[1].second;
  double O1 = parameters[2].first;
  double O2 = parameters[2].second;
  double X = point.first;
  double Y = point.second;
  return (std::pair<double, double>){ X1 * X + Y1 * Y + O1,
                                      X2 * X + Y2 * Y + O2 };
}

std::pair<double, double>
iterate (std::vector<transform> transformations,
         std::pair<double, double> point)
{
  return std::pair<double, double> ({ 0., 0. });
}