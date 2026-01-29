#include <cuda_runtime.h>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>

// CUDA kernel - now uses flat array directly
__global__ void evoDistanceKernel(
    const uint8_t* leafGenotype,
    const uint8_t* branchAllelesFlat,  // CHANGED: flat array instead of pointers
    int* distances,
    size_t numBranches,
    size_t L
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(idx >= numBranches) return;
    
    // Calculate offset into flat array
    const uint8_t* alleles = branchAllelesFlat + (idx * L);
    
    int dist = 0;
    
    #pragma unroll 4
    for(size_t i = 0; i < L; i++) {
        uint8_t leaf = leafGenotype[i];
        uint8_t branch = alleles[i];
        
        uint8_t first = branch & 0x0F;
        uint8_t second = branch >> 4;
        
        dist += (leaf != 0) & (leaf != first) & (leaf != second);
    }
    
    distances[idx] = dist;
}

// Updated host function - simpler signature
extern "C" void computeDistancesGPU(
    const uint8_t* d_leafGenotype,
    const uint8_t* d_branchAllelesFlat,  // CHANGED: flat array
    int* d_distances,
    size_t numBranches,
    size_t L
) {
    if(numBranches == 0) return;
    
    int threadsPerBlock = 256;
    int blocksPerGrid = (numBranches + threadsPerBlock - 1) / threadsPerBlock;
    
    evoDistanceKernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_leafGenotype,
        d_branchAllelesFlat,  // CHANGED
        d_distances,
        numBranches,
        L
    );
    
    cudaError_t launchErr = cudaGetLastError();
    if(launchErr != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel launch failed: ") + cudaGetErrorString(launchErr));
    }
    
    cudaError_t syncErr = cudaDeviceSynchronize();
    if(syncErr != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel execution failed: ") + cudaGetErrorString(syncErr));
    }
}