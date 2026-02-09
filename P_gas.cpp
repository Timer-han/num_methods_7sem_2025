#include "P_gas.h"

double rho(double x, double t)
{
  return exp (t) * (cos (3 * x * M_PI) + 1.5);
}

double u (double x, double t)
{
  return cos (2 * M_PI * t) * sin (4 * M_PI * x);
}

double rho_0 (double x)
{
  return cos(3 * M_PI * x) + 1.5;
}

double u_0 (double x)
{
  return sin(4 * M_PI * x);
}

double d_rho (double x, double t)
{
  return -3 * M_PI * exp (t) * sin (3 * x * M_PI);
}

double d_u_t (double x, double t)
{
  return  -2 * M_PI * sin (2 * M_PI * t) * sin (4 * M_PI * x);
}

double d_u_x (double x, double t)
{
  return 4 * M_PI * cos (2 * M_PI * t) * cos (4 * M_PI * x);
}

double d_u_xx (double x, double t)
{
  return -16 * M_PI * M_PI * cos(2 * M_PI * t) * sin(4 * M_PI * x);
}

double d_rho_t (double x, double t)
{
  return exp (t) * (cos (3 * x * M_PI) + 1.5);
}

double P_gas::get_Fn_m (double x, double t)
{
  double result = u (x, t) * d_u_x (x, t) + d_u_t (x, t) +
                  get_P(rho(x,t)) * (-3 * M_PI*sin(3 * M_PI*x)) / (cos(3*M_PI*x) + 1.5)
                  - mu * d_u_xx(x,t) / rho(x, t);

  return result;
}

double P_gas::get_F0 (double x, double t)
{
  return 1 + u (x,t) * (-3 * M_PI*sin(3 * M_PI*x)) / (cos(3*M_PI*x) + 1.5) + d_u_x(x,t);
}

double P_gas::get_P (double Gn_m)
{
  return mode ? p_ro : p_gamma * pow(Gn_m, p_gamma - 1);
}
