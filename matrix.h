#pragma once
#include "P_scheme.h"
#include "P_gas.h"
#include <vector>
#include <iostream>

enum class solver_mode : bool {V = true, G = false};
enum class debug_mode : bool {debug = true, no_debug = false};

class Matrix
{
  public :
  Matrix (P_gas gas, P_scheme scheme);

  P_gas gas;
  P_scheme scheme;

  int Dim;
  int step = 0;
  std::vector<double> matrix_vec_a;
  std::vector<double> matrix_vec_b;
  std::vector<double> matrix_vec_c;
  std::vector<double> rhs_vector;
  std::vector<double> solution_G;
  std::vector<double> solution_V;

  double get_mu ();
  double calc_nev_C (solver_mode mode);
  double calc_nev_l2 (solver_mode mode);
  double calc_nev_w2 (solver_mode mode, double res_l2);
  void calc_and_print_all_res(solver_mode mode);
  void calc_and_print_all_res_tex(solver_mode mode, char *filename);

  void init_matrix_G();
  void init_matrix_V();
  void solve (solver_mode mode, debug_mode debug = debug_mode::no_debug);

  void run_task ();
};

// Функции для вычисления норм
double norma_c(std::vector<double> vec, int Dim);
double norma_l(std::vector<double> vec, int Dim, double h);
double norma_w(std::vector<double> vec, int Dim, double h);

// Функции для вычисления невязок
double get_residual_C(std::vector<double> original, std::vector<double> divided, int orig_steps, int mult);
double get_residual_l2(std::vector<double> original, std::vector<double> divided, int orig_steps, int mult);
double get_residual_w2(std::vector<double> original, std::vector<double> divided, int orig_steps, int mult);
