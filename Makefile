# Makefile with platform-aware OpenMP enforcement and shared library output

UNAME_S := $(shell uname -s)

# Shared library extension per platform
ifeq ($(UNAME_S),Darwin)
  LIB_EXT := .dylib
  SHARED_FLAG := -dynamiclib
else ifeq ($(UNAME_S),Linux)
  LIB_EXT := .so
  SHARED_FLAG := -shared
else
  # Windows (MSYS2 / MinGW)
  LIB_EXT := .dll
  SHARED_FLAG := -shared
endif

# Compiler and flags per platform
ifeq ($(UNAME_S),Darwin)
  LLVM_DIR := $(shell brew --prefix llvm 2>/dev/null)
  ifneq ($(wildcard $(LLVM_DIR)/bin/clang++),)
    CXX := $(LLVM_DIR)/bin/clang++
    CXXFLAGS := -O3 -march=native -std=c++17 -fPIC -fopenmp \
                -I$(LLVM_DIR)/include -L$(LLVM_DIR)/lib
    LDFLAGS := -lomp -lpthread
    # Set install_name so the dylib can be found at load time from its own directory
    EXTRA_LDFLAGS := -install_name @rpath/liblinefit$(LIB_EXT)
  else
    $(error LLVM with OpenMP support not found. \
Install it using Homebrew:\n\n  brew install llvm\n\
Then run:\n  make\n\
You may also need to run:\n  export PATH="`brew --prefix llvm`/bin:$$PATH"\n)
  endif
else
  CXX := g++
  CXXFLAGS := -O3 -march=native -std=c++17 -fPIC -fopenmp
  LDFLAGS := -lpthread
  EXTRA_LDFLAGS :=
endif

TARGET := liblinefit$(LIB_EXT)
SRC := linefit.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SHARED_FLAG) $(SRC) -o $(TARGET) $(LDFLAGS) $(EXTRA_LDFLAGS)

clean:
	rm -f liblinefit.so liblinefit.dylib liblinefit.dll