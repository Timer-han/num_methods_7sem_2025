#pragma once
#include <cstdio>
#include "matrix.h"

inline void init_task_log(double mu, double C_rho, int mode, char* filename, solver_mode smode)
{
    FILE *file = fopen(filename, "w");
    if (!file) return;
    
    fprintf(file, "\\begin{table}[H]\n");
    fprintf(file, "\\centering\n");
    fprintf(file, "\\caption{Результаты для mu=%.3lf, C\\_rho=%.1lf, mode=%d, solver=%s}\n", 
            mu, C_rho, mode, (smode == solver_mode::G) ? "G" : "V");
    fprintf(file, "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}\n");
    fprintf(file, "\\hline\n");
    fprintf(file, "$\\tau$ & \\multicolumn{3}{c|}{$M_x=10$} & \\multicolumn{3}{c|}{$M_x=100$} & \\multicolumn{3}{c|}{$M_x=1000$} & \\multicolumn{3}{c|}{$M_x=10000$} \\\\ \\hline\n");
    fprintf(file, " & $C$ & $L_2$ & $W_2$ & $C$ & $L_2$ & $W_2$ & $C$ & $L_2$ & $W_2$ & $C$ & $L_2$ & $W_2$ \\\\ \\hline\n");
    
    fclose(file);
}

inline void end_task_log(char* filename)
{
    FILE *file = fopen(filename, "a");
    if (!file) return;
    
    fprintf(file, "\n\\end{tabular}\n");
    fprintf(file, "\\end{table}\n");
    
    fclose(file);
}
