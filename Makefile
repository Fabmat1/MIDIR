# Makefile with platform-aware OpenMP enforcement

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  LLVM_DIR := $(shell brew --prefix llvm 2>/dev/null)
  ifneq ($(wildcard $(LLVM_DIR)/bin/clang++),)
    CXX := $(LLVM_DIR)/bin/clang++
    CXXFLAGS := -O3 -march=native -std=c++17 -fopenmp -I$(LLVM_DIR)/include -L$(LLVM_DIR)/lib
    LDFLAGS := -lomp
  else
    $(error LLVM with OpenMP support not found. \
Install it using Homebrew:\n\n  brew install llvm\n\
Then run:\n  make\n\
You may also need to run:\n  export PATH="`brew --prefix llvm`/bin:$$PATH"\n)
  endif
else
  CXX := g++
  CXXFLAGS := -O3 -fopenmp -march=native -std=c++17
endif

TARGET = linefit
SRC = linefit.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
