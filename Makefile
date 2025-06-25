# Makefile
CXX = g++
CXXFLAGS = -O3 -fopenmp -march=native -std=c++17
TARGET = linefit
SRC = linefit.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TARGET)
