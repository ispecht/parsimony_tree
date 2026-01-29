#include <cuda_runtime.h>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>

// CUDA kernel for parallel distance computation
__global__ void evoDistanceKernel(
    const uint8_t* leafGenotype,
    uint8_t** branchAlleles,
    int* distances,
    size_t numBranches,
    size_t L
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(idx >= numBranches) return;
    
    const uint8_t* alleles = branchAlleles[idx];
    int dist = 0;
    
    for(size_t i = 0; i < L; i++) {
        uint8_t leaf = leafGenotype[i];
        uint8_t branch = alleles[i];
        
        uint8_t first = branch & 0x0F;
        uint8_t second = branch >> 4;
        
        dist += (leaf != 0) & (leaf != first) & (leaf != second);
    }
    
    distances[idx] = dist;
}

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
    
    // Check for kernel launch errors
    cudaError_t launchErr = cudaGetLastError();
    if(launchErr != cudaSuccess) {
        std::cerr << "Kernel launch error: " << cudaGetErrorString(launchErr) << std::endl;
        throw std::runtime_error(std::string("CUDA kernel launch failed: ") + cudaGetErrorString(launchErr));
    }
    
    // Wait for kernel to complete
    cudaError_t syncErr = cudaDeviceSynchronize();
    if(syncErr != cudaSuccess) {
        std::cerr << "Kernel sync error: " << cudaGetErrorString(syncErr) << std::endl;
        throw std::runtime_error(std::string("CUDA kernel execution failed: ") + cudaGetErrorString(syncErr));
    }
    
}