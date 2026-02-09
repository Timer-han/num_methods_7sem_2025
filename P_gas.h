#pragma once
#include <cmath>

class P_gas
{
  public :
  // Параметры сетки : T и Х
  double Segm_T;
  double Segm_X;
  // Параметры газа
  double p_ro; // С_rho для зависимости C_rho*rho
  double p_gamma; // gamma для rho^gamma
  double mu; // вязкость газа mu
  bool mode; // 0 - rho^gamma; 1 - C * rho

  P_gas (double Segm_T, double Segm_X, double p_ro, double p_gamma, double mu, bool mode) :
    Segm_T(Segm_T),
    Segm_X(Segm_X),
    p_ro(p_ro),
    p_gamma(p_gamma),
    mu(mu),
    mode(mode)
  {}

    double get_Fn_m (double x, double t);
    double get_F0 (double x, double t);
    double get_P (double Gn_m);
};

double rho (double x, double t);
double u (double x, double t);
double rho_0 (double x);
double u_0 (double x);
double d_rho (double x, double t);
double d_u_t (double x, double t);
double d_u_x (double x, double t);
double d_u_xx (double x, double t);
