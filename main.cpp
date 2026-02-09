#include <iostream>
#include "matrix.h"
#include "tex.h"
#include <sys/stat.h>

// Функция для создания директории, если её нет
void create_directory(const char* path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0700);
    }
}

int main(int argc, char const *argv[])
{
  // Создаём директорию для логов
  create_directory("./tex_logs");

  std::cout << "Запуск теста с упрощенными параметрами...\n" << std::endl;

  // Упрощенный тест с меньшим количеством итераций
  std::vector<int> mesh_steps = {10, 100};
  std::vector<int> modes = {0};
  std::vector<double> possible_C_rho = {1};
  std::vector<double> possible_mu = {0.1};

  for (int mode : modes)
    {
      for (double C_rho : possible_C_rho)
        {
          for (double mu : possible_mu)
            {
              char filename_G[100];
              char filename_V[100];
              sprintf(filename_G, "./tex_logs/mu_%.3lf_C_rho_%.1lf_mode%d_G.tex", mu, C_rho, mode);
              sprintf(filename_V, "./tex_logs/mu_%.3lf_C_rho_%.1lf_mode%d_V.tex", mu, C_rho, mode);
              
              init_task_log(mu, C_rho, mode, filename_G, solver_mode::G);
              init_task_log(mu, C_rho, mode, filename_V, solver_mode::V);

              for (int N : mesh_steps )
                {
                  FILE *file = fopen(filename_G, "a");
                  fprintf(file,"\n");
                  fprintf(file,"$%.4lf$ & ",1./N);
                  fclose(file);

                  file = fopen(filename_V, "a");
                  fprintf(file,"\n");
                  fprintf(file,"$%.4lf$ & ",1./N);
                  fclose(file);

                  for (int M_x : mesh_steps)
                    {
                      double tau = 1. / N;
                      double h_x = 1. / M_x;
                      P_gas gas (1., 1., C_rho, 1.4, mu, mode);
                      P_scheme scheme (M_x, N, h_x, tau);
                      Matrix matrix (gas, scheme);
                      printf("TASK: mode=%d C_rho=%.1lf, mu=%.3lf M_x=%d, N=%d\n", 
                             mode, C_rho, mu, M_x, N);
                      matrix.run_task();
                      matrix.calc_and_print_all_res_tex(solver_mode::G, filename_G);
                      matrix.calc_and_print_all_res_tex(solver_mode::V, filename_V);
                      printf("OK\n");
                    }

                    file = fopen(filename_G, "a");
                    fprintf(file," \\\\ \\hline");
                    fclose(file);

                    file = fopen(filename_V, "a");
                    fprintf(file," \\\\ \\hline");
                    fclose(file);
                }
              end_task_log(filename_G);
              end_task_log(filename_V);
            }
        }
      }

  // Тест на сходимость
  std::cout << "\n\nТест на сходимость с измельчением сетки...\n" << std::endl;
  
  for (int steps = 10; steps <= 100; steps*=10)
  {
    printf("\nSTEPS = %d\n", steps);
    
    // Запуск на стандартной сетке
    double tau = 1./steps;
    double h_x = 1./steps;

    P_gas gas1 (1., 1., 1, 1.4, 0.1, 0);
    P_scheme scheme1 (steps, steps, h_x, tau);
    Matrix matrix1(gas1, scheme1);
    matrix1.run_task();
    auto vec_g1 = matrix1.solution_G;
    auto vec_v1 = matrix1.solution_V;

    // Запуск с сеткой h/2
    tau /= 2;
    h_x /= 2;

    P_gas gas2 (1., 1., 1, 1.4, 0.1, 0);
    P_scheme scheme2 (steps * 2, steps * 2, h_x, tau);
    Matrix matrix2(gas2, scheme2);
    matrix2.run_task();
    auto vec_g2 = matrix2.solution_G;
    auto vec_v2 = matrix2.solution_V;

    // Запуск с сеткой h/4
    tau /= 2;
    h_x /= 2;

    P_gas gas3 (1., 1., 1, 1.4, 0.1, 0);
    P_scheme scheme3 (steps * 4, steps * 4, h_x, tau);
    Matrix matrix3(gas3, scheme3);
    matrix3.run_task();
    auto vec_g3 = matrix3.solution_G;
    auto vec_v3 = matrix3.solution_V;

    double res_v1_c = get_residual_C(vec_v1, vec_v2, steps, 2);
    double res_v1_l = get_residual_l2(vec_v1, vec_v2, steps, 2);
    double res_v1_w = get_residual_w2(vec_v1, vec_v2, steps, 2);

    double res_g1_c = get_residual_C(vec_g1, vec_g2, steps, 2);
    double res_g1_l = get_residual_l2(vec_g1, vec_g2, steps, 2);
    double res_g1_w = get_residual_w2(vec_g1, vec_g2, steps, 2);

    printf("RES_V1_C = %.3le, RES_V1_L2 = %.3le RES_V1_W2 = %.3le\n", res_v1_c, res_v1_l, res_v1_w);
    printf("RES_G1_C = %.3le, RES_G1_L2 = %.3le RES_G1_W2 = %.3le\n\n", res_g1_c, res_g1_l, res_g1_w);

    // Пересоздаём векторы для второго сравнения
    vec_g1 = matrix1.solution_G;
    vec_v1 = matrix1.solution_V;

    double res_v2_c = get_residual_C(vec_v1, vec_v3, steps, 4);
    double res_v2_l = get_residual_l2(vec_v1, vec_v3, steps, 4);
    double res_v2_w = get_residual_w2(vec_v1, vec_v3, steps, 4);

    double res_g2_c = get_residual_C(vec_g1, vec_g3, steps, 4);
    double res_g2_l = get_residual_l2(vec_g1, vec_g3, steps, 4);
    double res_g2_w = get_residual_w2(vec_g1, vec_g3, steps, 4);

    printf("RES_V2_C = %.3le, RES_V2_L2 = %.3le RES_V2_W2 = %.3le\n", res_v2_c, res_v2_l, res_v2_w);
    printf("RES_G2_C = %.3le, RES_G2_L2 = %.3le RES_G2_W2 = %.3le\n\n", res_g2_c, res_g2_l, res_g2_w);

    printf("ORIGINAL RES :\n");
    matrix1.calc_and_print_all_res(solver_mode::G);
    matrix1.calc_and_print_all_res(solver_mode::V);
    printf("-------------------------\n");
  }

  std::cout << "\n\nРабота программы завершена успешно!" << std::endl;
  std::cout << "Результаты сохранены в директории ./tex_logs/" << std::endl;

  return 0;
}
