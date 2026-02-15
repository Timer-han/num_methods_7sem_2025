#include "matrix.h"
#include <cmath>
#include <algorithm>

double find_norm_c(const std::vector<double> &v)
{
    double mx = 0;
    for (size_t i = 0; i < v.size(); i++)
        mx = std::max(mx, std::abs(v[i]));
    return mx;
}

double find_norm(const std::vector<double> &v)
{
    double mx = 0, h = 1. / (v.size() - 1);
    for (size_t i = 0; i < v.size(); i++)
        mx += v[i] * v[i];
    return std::sqrt(mx * h);
}

double find_norm_L2h(const std::vector<double> &v)
{
    double mx = find_norm(v), h = 1. / (v.size() - 1);
    mx = mx * mx / h;
    for (size_t i = 0; i < v.size(); i++)
        mx += 0.5 * v[i] * v[i];
    return std::sqrt(mx * h);
}

double find_norm_12(const std::vector<double> &v)
{
    double mx = 0, h = 1. / (v.size() - 1);
    for (size_t i = 0; i < v.size() - 1; i++)
        mx += v[i] * v[i];
    mx *= h;
    double L2h = find_norm_L2h(v);
    return std::sqrt(mx + L2h * L2h);
}