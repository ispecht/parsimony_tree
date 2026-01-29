# Compiler and flags
CXX := g++
NVCC := nvcc
CXXFLAGS := -std=c++20 -Wall -Wextra -O3 -march=native
NVCCFLAGS := -O3 -std=c++17 --compiler-options -fPIC

# Directories
SRC_DIR := src
INC_DIR := include
BUILD_DIR := build
BIN_DIR := bin

# Target executable
TARGET := $(BIN_DIR)/parsimony_tree

# Source files
CXX_SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
CU_SOURCES := $(wildcard $(SRC_DIR)/*.cu)
CXX_OBJECTS := $(CXX_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
CU_OBJECTS := $(CU_SOURCES:$(SRC_DIR)/%.cu=$(BUILD_DIR)/%.o)

OBJECTS := $(CXX_OBJECTS) $(CU_OBJECTS)

# CUDA paths (adjust if needed)
CUDA_PATH := /usr/local/cuda
CUDA_INC := $(CUDA_PATH)/include
CUDA_LIB := $(CUDA_PATH)/lib64
CUDA_LIBS := -lcudart

# Default target
all: directories $(TARGET)

# Create necessary directories
directories:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Link executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ -L$(CUDA_LIB) $(CUDA_LIBS)
	@echo "Built: $(TARGET)"

# Compile C++ source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -I$(CUDA_INC) -c $< -o $@

# Compile CUDA source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -I$(INC_DIR) -c $< -o $@

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean directories