# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O3 -march=native

# CUDA compiler - try to find nvcc automatically
NVCC := $(shell which nvcc 2>/dev/null || echo /usr/local/cuda/bin/nvcc)
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

# CUDA paths - auto-detect from nvcc location
CUDA_PATH := $(shell dirname $(shell dirname $(NVCC)) 2>/dev/null || echo /usr/local/cuda)
CUDA_INC := $(CUDA_PATH)/include
CUDA_LIB := $(CUDA_PATH)/lib64
CUDA_LIBS := -lcudart

# Dependency files
DEPS := $(OBJECTS:.o=.d)

# Include directories for all compilation
INCLUDES := -I$(INC_DIR) -I$(CUDA_INC)

# Default target
all: directories $(TARGET)

# Create necessary directories
directories:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Link executable
$(TARGET): $(OBJECTS)
	@echo "Linking $@..."
	@$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ -L$(CUDA_LIB) $(CUDA_LIBS)
	@echo "Built: $(TARGET)"

# Compile C++ source files with dependency generation
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

# Compile CUDA source files with dependency generation
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	@echo "Compiling CUDA $<..."
	@$(NVCC) $(NVCCFLAGS) $(INCLUDES) -M $< -o $(@:.o=.d) -odir $(BUILD_DIR)
	@$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@

# Include dependency files
-include $(DEPS)

# Clean build artifacts
clean:
	@echo "Cleaning..."
	@rm -rf $(BUILD_DIR) $(BIN_DIR)
	@echo "Clean complete"

# Clean and rebuild
rebuild: clean all

# Print configuration (useful for debugging)
info:
	@echo "Configuration:"
	@echo "  CXX:        $(CXX)"
	@echo "  NVCC:       $(NVCC)"
	@echo "  CUDA_PATH:  $(CUDA_PATH)"
	@echo "  CUDA_INC:   $(CUDA_INC)"
	@echo "  CUDA_LIB:   $(CUDA_LIB)"
	@echo "  CXX Sources: $(CXX_SOURCES)"
	@echo "  CU Sources:  $(CU_SOURCES)"
	@echo "  Objects:     $(OBJECTS)"

# Phony targets
.PHONY: all clean directories rebuild info