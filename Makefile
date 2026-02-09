CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
TARGET = evm
TEST_TARGET = evm_test
FULL_TARGET = evm_full
SOURCES = main.cpp matrix.cpp P_gas.cpp Residual.cpp
TEST_SOURCES = main_test.cpp matrix.cpp P_gas.cpp Residual.cpp
FULL_SOURCES = main_full_latex.cpp matrix.cpp P_gas.cpp Residual.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) -lm

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Быстрый тест генерации таблиц (2 параметра mu, tau, h)
test:
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_SOURCES) -lm
	./$(TEST_TARGET)

# Полная генерация таблиц (все параметры, занимает несколько часов)
full:
	$(CXX) $(CXXFLAGS) -o $(FULL_TARGET) $(FULL_SOURCES) -lm
	./$(FULL_TARGET)

# Запуск полной генерации в фоновом режиме
background:
	$(CXX) $(CXXFLAGS) -o $(FULL_TARGET) $(FULL_SOURCES) -lm
	nohup ./$(FULL_TARGET) > full_log.txt 2>&1 &
	@echo "Процесс запущен в фоне. Проверяйте прогресс: tail -f full_log.txt"

clean:
	rm -f $(OBJECTS) $(TARGET) $(TEST_TARGET) $(FULL_TARGET)
	rm -rf tex_logs results full_log.txt

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run test full background
