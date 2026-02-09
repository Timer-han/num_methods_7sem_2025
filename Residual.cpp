#include "matrix.h"
#include <cmath>

double norma_c(std::vector<double> vec, int Dim)
{
  int m;
  double n_c, tmp;
  n_c = fabs (vec[0]);
  for(m=1; m<Dim; m++)
     {
        tmp = fabs (vec[m]);
        if (n_c < tmp) n_c=tmp;
     }
  return n_c;
}

double norma_l(std::vector<double> vec,int Dim, double h)
{
  int m;
  double n_l, tmp;
  n_l = vec[0] * vec[0];
  for (m = 1; m < Dim; m++)
     {
        n_l += vec[m] * vec[m];
     }
  return sqrt (h * n_l);
}

double norma_w(std::vector<double> vec,int Dim,double h)
{
  int m;
  double n_l, n_w, tmp;
  n_l = vec[0]*vec[0];
  n_w = 0;
  for (m = 1; m < Dim; m++)
     {
        n_l += vec[m] * vec[m];
        n_w += (vec[m] - vec[m - 1]) * (vec[m] - vec[m - 1]);
     }
  return sqrt (h * n_l + n_w / h);
}

double get_residual_C(std::vector<double> original, std::vector<double> divided, int orig_steps, int mult)
{
  for (int i = 0; i <= orig_steps; i++)
    {
      original[i] -= divided[i * mult];
    }

    return norma_c(original, orig_steps+1);
}

double get_residual_l2(std::vector<double> original, std::vector<double> divided, int orig_steps, int mult)
{
  for (int i = 0; i <= orig_steps; i++)
    {
      original[i] -= divided[i * mult];
    }

    return norma_l(original, orig_steps+1, 1./orig_steps);
}

double get_residual_w2(std::vector<double> original, std::vector<double> divided, int orig_steps, int mult)
{
  for (int i = 0; i <= orig_steps; i++)
    {
      original[i] -= divided[i * mult];
    }

    return norma_w(original, orig_steps+1, 1./orig_steps);
}
