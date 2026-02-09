#include "matrix.h"
#include <cstdio>

Matrix::Matrix (P_gas gas, P_scheme scheme) :
gas(gas), scheme(scheme)
{
  Dim = scheme.Dim;
  matrix_vec_a.resize(Dim);
  matrix_vec_b.resize(Dim);
  matrix_vec_c.resize(Dim);
  rhs_vector.resize(Dim);
  solution_G.resize(Dim);
  solution_V.resize(Dim);
  for (int i = 0; i < Dim; i++)
    {
      double x = scheme.get_point_x_by_i(i);
      solution_G[i] = log(rho_0(x));
      solution_V[i] = u_0(x);
    }
  step = 0;
}

double Matrix::get_mu ()
{
  double max = -1;

  for (int i = 0; i < Dim; i++)
    {
      double candidate = exp(-solution_G[i]);
      if (candidate > max)
      {
        max = candidate;
      }
    }
  return gas.mu * max;
}

void Matrix::init_matrix_G ()
{
  double tau = scheme.tau;
  double h_x = scheme.h_x;
  solution_V[0] = 0;
  solution_V[Dim-1] = 0;

  // основная часть
  for (int m = 1; m < Dim-1; m++)
    {
      double Vn_m = solution_V[m];
      double Vn_m_r = solution_V[m+1];
      double Vn_m_l = solution_V[m-1];
      matrix_vec_a[m] = 1;
      matrix_vec_b[m] = (tau / (4 * h_x)) * (Vn_m + Vn_m_r);
      matrix_vec_c[m] = -(tau / (4 * h_x)) * (Vn_m + Vn_m_l);

      double Gn_m = solution_G[m];
      rhs_vector[m] = Gn_m + (Gn_m - 2) * (Vn_m_r - Vn_m_l) * tau / (4 * h_x) + tau*gas.get_F0(h_x * m, step * tau);
    }

  //граничные условия : 0 строка
  double Vn_0 = solution_V[0];
  double Vn_1 = solution_V[1];
  double Vn_2 = solution_V[2];
  double Vn_3 = solution_V[3];
  double Gn_0 = solution_G[0];
  double Gn_1 = solution_G[1];
  double Gn_2 = solution_G[2];
  double Gn_3 = solution_G[3];

  matrix_vec_c[0] = 0;
  matrix_vec_a[0] = (1  - (Vn_0 * tau) / (2 * h_x));
  matrix_vec_b[0] = (Vn_1 * tau ) / (2 * h_x);
  rhs_vector[0] = Gn_0 - tau * (2 - Gn_0) * (Vn_1 - Vn_0) / (2 * h_x) +
                  (tau / (4 * h_x)) * (-Vn_3 * Gn_3 + 4 * Vn_2 * Gn_2 - 5 * Vn_1 * Gn_1 + 2 * Vn_0 * Gn_0 +
                  (2 - Gn_0) * (-Vn_3 + 4 * Vn_2 - 5 * Vn_1 + 2 * Vn_0)) + tau*gas.get_F0(0, step * tau);

  // граничные условия : строка Dim-1
  double Vn_M_0 = solution_V[Dim - 1];
  double Vn_M_1 = solution_V[Dim - 2];
  double Vn_M_2 = solution_V[Dim - 3];
  double Vn_M_3 = solution_V[Dim - 4];
  double Gn_M_0 = solution_G[Dim - 1];
  double Gn_M_1 = solution_G[Dim - 2];
  double Gn_M_2 = solution_G[Dim - 3];
  double Gn_M_3 = solution_G[Dim - 4];

  matrix_vec_b[Dim - 1] = 0;
  matrix_vec_a[Dim - 1] = (1 + (Vn_M_0 * tau) / (2 * h_x));
  matrix_vec_c[Dim - 1] = (-Vn_M_1 * tau) / (2 * h_x);
  rhs_vector[Dim - 1] = Gn_M_0 - (tau  / (2 * h_x)) * (2 - Gn_M_0) * (Vn_M_0 - Vn_M_1) -
                        (tau / (4 * h_x)) * (2 * Gn_M_0 * Vn_M_0 - 5 * Gn_M_1 * Vn_M_1 + 4 * Gn_M_2 * Vn_M_2 - Gn_M_3 * Vn_M_3 +
                        (2 - Gn_M_0) * (2 * Vn_M_0 - 5 * Vn_M_1 + 4 * Vn_M_2 - Vn_M_3 )) + tau*gas.get_F0(h_x * (Dim - 1), step * tau);
}

void Matrix::init_matrix_V ()
{
  // первая строка фиктивная
  matrix_vec_a[0] = 1;
  matrix_vec_b[0] = 0;
  matrix_vec_c[0] = 0;
  rhs_vector[0] = 0;

  // основная часть
  double tau = scheme.tau;
  double h_x = scheme.h_x;
  double mu = get_mu();

  for (int m = 1; m < Dim-1; m++)
    {
      auto Vn_m = solution_V[m];
      auto Vn_m_r = solution_V[m + 1];
      auto Vn_m_l = solution_V[m - 1];

      matrix_vec_a[m] = 1 + (2 * mu * tau) / (h_x * h_x);
      matrix_vec_b[m] = tau * (Vn_m + Vn_m_r) / (6 * h_x) - tau * mu / (h_x * h_x);
      matrix_vec_c[m] = -tau *(Vn_m + Vn_m_l) / (6 * h_x) - tau * mu / (h_x * h_x);

      auto Gn_m = solution_G[m];
      auto Gn_m_r = solution_G[m + 1];
      auto Gn_m_l = solution_G[m - 1];

      auto P_G = gas.get_P(exp(Gn_m));
      auto Fn_m = gas.get_Fn_m (h_x * m, tau * step);

      rhs_vector[m] = Vn_m + tau * Fn_m - tau * P_G * (Gn_m_r - Gn_m_l) / (2 * h_x)
                      - tau * (mu - gas.mu * exp(-Gn_m)) * (Vn_m_r - 2*Vn_m + Vn_m_l) / (h_x * h_x);
    }

  // последняя строка тоже фиктивная
  matrix_vec_a[Dim - 1] = 1;
  matrix_vec_b[Dim - 1] = 0;
  matrix_vec_c[Dim - 1] = 0;
  rhs_vector[Dim - 1] = 0;
}

void Matrix::solve (solver_mode mode, debug_mode debug)
{
  std::vector<double> & solution = (mode == solver_mode::V) ? solution_V : solution_G;
  // Метод прогонки : вперед

  double znam = matrix_vec_a[0];
  matrix_vec_a[0] =  (-1 * matrix_vec_b[0]) / znam;
  matrix_vec_b[0] = rhs_vector[0] / znam;

  for (int i = 1; i < Dim - 1; i++)
    {
      znam = matrix_vec_a[i] + matrix_vec_c[i] * matrix_vec_a[i-1];
      matrix_vec_a[i] = (-1 * matrix_vec_b[i]) / znam;
      double chis = rhs_vector[i] - matrix_vec_c[i] * matrix_vec_b[i-1];
      matrix_vec_b[i] = chis / znam;
    }
  // Метод прогонки : обратно
  solution[Dim - 1] = (rhs_vector[Dim - 1] - matrix_vec_c[Dim - 1]*matrix_vec_b[Dim - 2]) /
                      (matrix_vec_a[Dim - 1] + matrix_vec_c[Dim - 1] * matrix_vec_a[Dim - 2]);
  for (int i = Dim - 2; i >=0; i--)
    {
      solution[i] = matrix_vec_a[i] * solution[i+1] + matrix_vec_b[i];
    }

    if (debug == debug_mode::debug)
     {
       for (int i = 0; i < Dim; i++)
         {
           double x = i * scheme.h_x;
           double t = step * scheme.tau;

           if (mode == solver_mode::G) // Vector G
             solution[i] = log (rho (x, t));
           else
             solution[i] = u (x, t);
         }
     }
}

void Matrix::run_task ()
{
  step = 1;
  init_matrix_G ();
  solve (solver_mode::G);
  init_matrix_V ();
  solve (solver_mode::V);

  // step 1 passed go futher
  for (int i = 2; i <= scheme.N; i++)
    {
      step ++;
      init_matrix_G ();
      solve (solver_mode::G);
      init_matrix_V ();
      solve (solver_mode::V);
    }
}

double Matrix::calc_nev_C (solver_mode mode)
{
  double max = 0;
  int max_i = -1;
  for (int i = 0; i < Dim; i++)
    {
      double x = scheme.h_x * i;
      double tau = scheme.tau * step;
      if (mode == solver_mode::G)
        {
          if (fabs(solution_G[i] - log (rho (x,tau))) > max)
            {
              max = fabs(solution_G[i] - log (rho (x,tau)));
              max_i = i;
            }
        }
        else
          {
            if (fabs (solution_V[i] - u (x, tau)) > max)
              {
                max = fabs (solution_V[i] - u (x, tau));
                max_i = i;
              }
          }
    }

  return max;
}

double Matrix::calc_nev_l2 (solver_mode mode)
{
  double sum = 0;
  for (int i = 0; i < Dim; i++)
    {
      double x = scheme.h_x * i;
      double tau = scheme.tau * step;

      if (mode == solver_mode::G)
        sum += (solution_G[i] - log (rho (x,tau))) *
               (solution_G[i] - log (rho (x,tau)));
        else
          sum += (solution_V[i] - u (x, tau)) * (solution_V[i] - u (x, tau));
    }
    sum = sum * scheme.h_x;
    return sqrt(sum);
}

double Matrix::calc_nev_w2 (solver_mode mode, double res_l2)
{
  double sum = 0;
  double tau = scheme.tau * step;
  double prev_res = (mode == solver_mode::G) ? (solution_G[0] - log (rho (0, tau))) :
                                               (solution_V[0] - u (0, tau));
  for (int i = 1; i < Dim; i++)
    {
      double x = scheme.h_x * i;
      double res{};

      if (mode == solver_mode::G)
        res = (solution_G[i] - log (rho (x,tau)));
      else
        res = (solution_V[i] - u (x, tau));

      sum += (res - prev_res) *
             (res - prev_res);
      prev_res = res;
    }
    sum = res_l2 * res_l2 + sum / scheme.h_x;
    return sqrt (sum);
}

void Matrix::calc_and_print_all_res(solver_mode mode)
{
  double C_res = calc_nev_C (mode);
  double l2_res = calc_nev_l2 (mode);
  double w2_res = calc_nev_w2 (mode, l2_res);
  printf ("\nMODE %d step %d C_RES = %lf L2_RES = %lf W2_RES = %lf\n",
          static_cast<bool>(mode), step, C_res, l2_res, w2_res);
}

void Matrix::calc_and_print_all_res_tex(solver_mode mode, char *filename)
{
  double C_res = calc_nev_C (mode);
  double l2_res = calc_nev_l2 (mode);
  double w2_res = calc_nev_w2 (mode, l2_res);
  FILE *file = fopen(filename, "a");
  fprintf(file, "$%.3le$ $%.3le$ $%.3le$ ",C_res, l2_res, w2_res);
  if (Dim -1 != 10000)
    fprintf(file,"&");
  fclose(file);
}
