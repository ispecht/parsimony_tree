#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct InitBranch;

struct InitNode {
    std::vector<InitBranch*> childBranches;
    InitBranch* parentBranch;
    std::string name;
    std::string seq = "";
};

struct InitBranch {
    InitNode* parent;
    InitNode* child;
    std::vector<uint8_t> alleles;
    int nMuts;
};

// Update the function signature
extern "C" void computeDistancesGPU(
    const uint8_t* d_leafGenotype,
    const uint8_t* d_branchAllelesFlat,  // CHANGED: no longer **
    int* d_distances,
    size_t numBranches,
    size_t L
);