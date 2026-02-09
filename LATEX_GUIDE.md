# Генерация LaTeX таблиц для задачи 3.11

## Описание задачи

Ваша задача (3.11) использует следующие уравнения:

```
G_t + 0.5(V Ĝ_x̊ + (V Ĝ)_x̊ + (2 - G)V_x̊) = 0
V_t + (1/3)(V V̂_x̊ + (V V̂)_x̊) + p̃'(Ĝ)Ĝ_x̊ = μ̃ V̂_xx̄ - (μ̃ - μe^(-Ĝ))V_xx̄ + f
```

## Доступные программы

В проекте есть три версии программы для генерации таблиц:

### 1. main_test.cpp (БЫСТРЫЙ ТЕСТ)
**Параметры:**
- mu: {0.1}
- C_rho: {1, 10}
- mode: {1} (только p(ρ) = C·ρ)
- tau: {0.1, 0.01}
- h: {0.1, 0.01}

**Время работы:** ~5-10 секунд

**Использование:**
```bash
make test
# или
g++ -std=c++11 -O2 -o evm_test main_test.cpp matrix.cpp P_gas.cpp Residual.cpp -lm
./evm_test
```

**Результаты:** `results/test_tables_G.tex` и `results/test_tables_V.tex`

### 2. main_full_latex.cpp (ПОЛНЫЙ НАБОР)
**Параметры:**
- mu: {0.1, 0.01, 0.001}
- C_rho: {1, 10, 100}
- mode: {0, 1} (оба типа p(ρ))
- tau: {0.1, 0.01, 0.001, 0.0001}
- h: {0.1, 0.01, 0.001, 0.0001}

**Количество вычислений:** 3×3×2×4×4 = 288 комбинаций

**Время работы:** 2-6 часов (зависит от процессора)

**Использование:**
```bash
make full
# или
g++ -std=c++11 -O2 -o evm_full main_full_latex.cpp matrix.cpp P_gas.cpp Residual.cpp -lm
./evm_full > full_log.txt 2>&1 &
```

**Результаты:** `results/full_tables_G.tex` и `results/full_tables_V.tex`

### 3. Запуск в фоновом режиме (рекомендуется)

Для длительных вычислений используйте:

```bash
nohup ./evm_full > full_log.txt 2>&1 &
```

Отслеживание прогресса:
```bash
tail -f full_log.txt
```

## Структура выходных таблиц

Каждая таблица имеет формат:

```latex
\begin{tabular}{ |l|l|l|l|l| }
\hline
\multicolumn{5}{|c|}{$\mu = 0.1, p(\rho) = 1\rho$} \\
\hline
$\tau\setminus h$ & $0.1$ & $0.01$ & $0.001$ & $0.0001$\\
\hline
$0.1$ & <C1> & <C2> & <C3> & <C4> \\
& <L2_1> & <L2_2> & <L2_3> & <L2_4> \\
& <W2_1> & <W2_2> & <W2_3> & <W2_4> \\
& <time1> & <time2> & <time3> & <time4> \\
\hline
...
\end{tabular}
```

Где:
- **C** - невязка в равномерной норме (максимум модуля)
- **L2** - невязка в интегральной L2 норме
- **W2** - невязка в норме Соболева W2
- **time** - произведение τ × h

## Результаты для разных параметров

### mode = 0: p(ρ) = C_rho × ρ^1.4
Уравнение состояния политропного газа

### mode = 1: p(ρ) = C_rho × ρ
Линейная зависимость давления от плотности

## Makefile команды

Добавьте в Makefile:

```makefile
test:
	g++ -std=c++11 -O2 -o evm_test main_test.cpp matrix.cpp P_gas.cpp Residual.cpp -lm
	./evm_test

full:
	g++ -std=c++11 -O2 -o evm_full main_full_latex.cpp matrix.cpp P_gas.cpp Residual.cpp -lm
	./evm_full

background:
	g++ -std=c++11 -O2 -o evm_full main_full_latex.cpp matrix.cpp P_gas.cpp Residual.cpp -lm
	nohup ./evm_full > full_log.txt 2>&1 &
	@echo "Процесс запущен в фоне. Проверяйте прогресс: tail -f full_log.txt"
```

## Интеграция в LaTeX документ

Скопируйте содержимое файлов в ваш LaTeX документ:

```latex
\documentclass{article}
\usepackage{diagbox} % Для диагональных заголовков
\usepackage{booktabs}

\begin{document}

\section{Результаты для G (логарифм плотности)}
\input{results/full_tables_G.tex}

\section{Результаты для V (скорость)}
\input{results/full_tables_V.tex}

\end{document}
```

## Примечания

1. **Память:** Для больших сеток (h=0.0001, N=10000) требуется достаточно RAM

2. **Точность:** Результаты выводятся в экспоненциальной форме с 6 значащими цифрами

3. **Ошибки:** Если вычисление не сходится, в таблице будет значение -1.000000e+00

4. **Проверка результатов:** Невязки должны уменьшаться при уменьшении h и τ

## Пример работы

```bash
cd evm_project

# Быстрый тест
make test
# Результат в results/test_tables_G.tex и results/test_tables_V.tex

# Полный расчёт (запуск в фоне)
make background
# Или вручную:
# g++ -std=c++11 -O2 -o evm_full main_full_latex.cpp matrix.cpp P_gas.cpp Residual.cpp -lm
# nohup ./evm_full > full_log.txt 2>&1 &

# Проверка прогресса
tail -f full_log.txt

# После завершения
ls -lh results/
cat results/full_tables_G.tex
```

## Настройка параметров

Если нужен другой набор параметров, отредактируйте в main_full_latex.cpp:

```cpp
std::vector<double> possible_mu = {0.1, 0.01, 0.001};       // Вязкость
std::vector<double> possible_C_rho = {1, 10, 100};          // Коэффициент в p(ρ)
std::vector<int> modes = {0, 1};                             // 0: ρ^γ, 1: C·ρ
std::vector<double> tau_values = {0.1, 0.01, 0.001, 0.0001};// Шаги по времени
std::vector<double> h_values = {0.1, 0.01, 0.001, 0.0001};  // Шаги по пространству
```

## Анализ результатов

Ожидаемое поведение:
- Невязки уменьшаются при измельчении сетки
- При фиксированном h невязки могут сначала уменьшаться, затем стабилизироваться
- Для стабильности нужно τ ≤ Ch² (условие Куранта)
- Последняя строка каждой группы (τ×h) показывает произведение шагов
