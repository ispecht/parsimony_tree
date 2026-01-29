#include <cuda_runtime.h>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <string>

// Forward declaration for InitBranch (we only need the structure layout)
struct InitBranch {
    void* parent;
    void* child;
    // We'll access alleles directly via pointer
    int nMuts;
};

// CUDA kernel for parallel distance computation
__global__ void evoDistanceKernel(
    const uint8_t* leafGenotype,
    uint8_t** branchAlleles,     // Array of pointers to branch alleles
    int* distances,
    size_t numBranches,
    size_t L
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(idx >= numBranches) return;
    
    const uint8_t* alleles = branchAlleles[idx];
    int dist = 0;
    
    #pragma unroll 4
    for(size_t i = 0; i < L; i++) {
        uint8_t leaf = leafGenotype[i];
        uint8_t branch = alleles[i];
        
        uint8_t first = branch & 0x0F;
        uint8_t second = branch >> 4;
        
        // Branchless distance calculation
        dist += (leaf != 0) & (leaf != first) & (leaf != second);
    }
    
    distances[idx] = dist;
}

// Host function to launch kernel (extern "C" for C++ linkage)
extern "C" void computeDistancesGPU(
    const uint8_t* d_leafGenotype,
    uint8_t** d_branchAlleles,
    int* d_distances,
    size_t numBranches,
    size_t L
) {
    if(numBranches == 0) return;
    
    // Launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (numBranches + threadsPerBlock - 1) / threadsPerBlock;
    
    evoDistanceKernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_leafGenotype,
        d_branchAlleles,
        d_distances,
        numBranches,
        L
    );
    
    // Wait for kernel to complete
    cudaError_t err = cudaDeviceSynchronize();
    if(err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA error: ") + cudaGetErrorString(err));
    }
}