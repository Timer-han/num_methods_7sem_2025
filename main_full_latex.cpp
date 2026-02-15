#include <iostream>
#include <fstream>
#include "matrix.h"
#include <sys/stat.h>

void create_directory(const char* path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0700);
    }
}

void write_table_header(std::ofstream &file, double mu, double C_rho, int mode) {
    file << "\\begin{tabular}{ |l|l|l|l|l| }\n";
    file << "\\hline\n";
    
    if (mode == 0) {
        file << "\\multicolumn{5}{|c|}{$\\mu = " << mu << ", p(\\rho) = " << C_rho << "\\rho^{1.4}$} \\\\\n";
    } else {
        file << "\\multicolumn{5}{|c|}{$\\mu = " << mu << ", p(\\rho) = " << C_rho << "\\rho$} \\\\\n";
    }
    
    file << "\\hline\n";
    file << "$\\tau\\setminus h$ & $0.1$ & $0.01$ & $0.001$ & $0.0001$\\\\\n";
    file << "\\hline\n";
}

void write_table_footer(std::ofstream &file) {
    file << "\\end{tabular}\n\n\n";
}

int main(int argc, char const *argv[])
{
    create_directory("./results");
    
    std::ofstream result_file("./results/full_tables_G.tex");
    std::ofstream result_file_V("./results/full_tables_V.tex");
    
    std::cout << "=== Генерация полных таблиц для задачи 3.11 ===\n" << std::endl;
    std::cout << "ВНИМАНИЕ: Полный запуск может занять несколько часов!\n" << std::endl;
    
    // ПОЛНЫЕ параметры
    std::vector<double> possible_mu = {0.1, 0.01, 0.001};
    std::vector<double> possible_C_rho = {1, 10, 100};
    std::vector<int> modes = {0, 1}; // 0 - rho^gamma, 1 - C*rho
    std::vector<double> tau_values = {0.1, 0.01, 0.001, 0.0001};
    std::vector<double> h_values = {0.1, 0.01, 0.001, 0.0001};
    
    for (double mu : possible_mu) {
        for (double C_rho : possible_C_rho) {
            for (int mode : modes) {
                
                std::cout << "\n========================================" << std::endl;
                std::cout << "mu = " << mu << ", C_rho = " << C_rho << ", mode = " << mode << std::endl;
                std::cout << "========================================\n" << std::endl;
                
                write_table_header(result_file, mu, C_rho, mode);
                write_table_header(result_file_V, mu, C_rho, mode);
                
                for (double tau : tau_values) {
                    int N = (int)(1.0 / tau);
                    
                    // Собираем все результаты для данного tau
                    std::vector<double> C_res_G_vec, l2_res_G_vec, l2h_res_G_vec, w2_res_G_vec;
                    std::vector<double> C_res_V_vec, l2_res_V_vec, l2h_res_V_vec, w2_res_V_vec;

                    for (double h_x : h_values) {
                        int M_x = (int)(1.0 / h_x);
                        
                        std::cout << "  tau = " << tau << ", h = " << h_x 
                                  << " (N=" << N << ", M_x=" << M_x << ")" << std::flush;
                        
                        try {
                            P_gas gas(1., 1., C_rho, 1.4, mu, mode);
                            P_scheme scheme(M_x, N, h_x, tau);
                            Matrix matrix(gas, scheme);
                            
                            matrix.run_task();
                            
                            double C_res_G = matrix.calc_nev_C(solver_mode::G);
                            double l2_res_G = matrix.calc_nev_l2(solver_mode::G);
                            double l2h_res_G = matrix.calc_nev_l2h(solver_mode::G);
                            double w2_res_G = matrix.calc_nev_w2(solver_mode::G);

                            double C_res_V = matrix.calc_nev_C(solver_mode::V);
                            double l2_res_V = matrix.calc_nev_l2(solver_mode::V);
                            double l2h_res_V = matrix.calc_nev_l2h(solver_mode::V);
                            double w2_res_V = matrix.calc_nev_w2(solver_mode::V);

                            C_res_G_vec.push_back(C_res_G);
                            l2_res_G_vec.push_back(l2_res_G);
                            l2h_res_G_vec.push_back(l2h_res_G);
                            w2_res_G_vec.push_back(w2_res_G);

                            C_res_V_vec.push_back(C_res_V);
                            l2_res_V_vec.push_back(l2_res_V);
                            l2h_res_V_vec.push_back(l2h_res_V);
                            w2_res_V_vec.push_back(w2_res_V);

                            std::cout << " ✓" << std::endl;
                            
                        } catch (...) {
                                std::cout << " ✗ (error)" << std::endl;
                                C_res_G_vec.push_back(-1);
                                l2_res_G_vec.push_back(-1);
                                l2h_res_G_vec.push_back(-1);
                                w2_res_G_vec.push_back(-1);
                                C_res_V_vec.push_back(-1);
                                l2_res_V_vec.push_back(-1);
                                l2h_res_V_vec.push_back(-1);
                                w2_res_V_vec.push_back(-1);
                        }
                    }
                    
                    // Записываем 4 строки: tau & C1 & C2 & ... \\
                    result_file.precision(6);
                    result_file_V.precision(6);
                    
                    // Строка 1: tau и C-нормы
                    // Строка 1: C-нормы
                    result_file << "$" << tau << "$ ";
                    result_file_V << "$" << tau << "$ ";
                    for (size_t i = 0; i < C_res_G_vec.size(); i++)
                    {
                        result_file << "& $" << std::scientific << C_res_G_vec[i] << "$ ";
                        result_file_V << "& $" << std::scientific << C_res_V_vec[i] << "$ ";
                    }
                    result_file << "\\\\\n";
                    result_file_V << "\\\\\n";

                    // Строка 2: L2-нормы
                    result_file << "& ";
                    result_file_V << "& ";
                    for (size_t i = 0; i < l2_res_G_vec.size(); i++)
                    {
                        result_file << "$" << std::scientific << l2_res_G_vec[i] << "$ ";
                        result_file_V << "$" << std::scientific << l2_res_V_vec[i] << "$ ";
                        if (i < l2_res_G_vec.size() - 1)
                        {
                            result_file << "& ";
                            result_file_V << "& ";
                        }
                    }
                    result_file << "\\\\\n";
                    result_file_V << "\\\\\n";

                    // Строка 3: L2h-нормы
                    result_file << "& ";
                    result_file_V << "& ";
                    for (size_t i = 0; i < l2h_res_G_vec.size(); i++)
                    {
                        result_file << "$" << std::scientific << l2h_res_G_vec[i] << "$ ";
                        result_file_V << "$" << std::scientific << l2h_res_V_vec[i] << "$ ";
                        if (i < l2h_res_G_vec.size() - 1)
                        {
                            result_file << "& ";
                            result_file_V << "& ";
                        }
                    }
                    result_file << "\\\\\n";
                    result_file_V << "\\\\\n";

                    // Строка 4: W2-нормы
                    result_file << "& ";
                    result_file_V << "& ";
                    for (size_t i = 0; i < w2_res_G_vec.size(); i++)
                    {
                        result_file << "$" << std::scientific << w2_res_G_vec[i] << "$ ";
                        result_file_V << "$" << std::scientific << w2_res_V_vec[i] << "$ ";
                        if (i < w2_res_G_vec.size() - 1)
                        {
                            result_file << "& ";
                            result_file_V << "& ";
                        }
                    }
                    result_file << "\\\\\n";
                    result_file_V << "\\\\\n";

                    result_file << "\\hline\n";
                    result_file_V << "\\hline\n";
                }
                
                write_table_footer(result_file);
                write_table_footer(result_file_V);
            }
        }
    }
    
    result_file.close();
    result_file_V.close();
    
    std::cout << "\n=== Генерация завершена ===\n" << std::endl;
    std::cout << "Результаты сохранены в:" << std::endl;
    std::cout << "  - ./results/full_tables_G.tex (для G)" << std::endl;
    std::cout << "  - ./results/full_tables_V.tex (для V)" << std::endl;
    
    return 0;
}
