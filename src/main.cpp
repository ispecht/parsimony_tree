#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <chrono>
#include <cuda_runtime.h>
#include "helpers.h"
#include "debug.h"
#include "structs.h"


// Compute min distance from a fixed genotype (leaf) to all possible genotypes allowed on a branch
int evoDistance(
    const std::vector<uint8_t>& leafGenotype,
    const std::vector<uint8_t>& alleles
) {
    int dist = 0;

    for(size_t i = 0; i < leafGenotype.size(); i++) {
        uint8_t leaf = leafGenotype[i];
        uint8_t branch = alleles[i];
        
        // Branchless: dist += (leaf != 0) && (leaf != low_nibble) && (leaf != high_nibble)
        uint8_t first = branch & 0x0F;
        uint8_t second = branch >> 4;
        
        dist += (leaf != 0) & (leaf != first) & (leaf != second);
    }

    return dist;
}

// Split branch, and return the new node in the middle
InitNode* splitBranch(
    InitBranch* branch,
    const std::vector<bool>& isLeft,
    std::vector<InitBranch*>& allBranches
) {

    // Rewiring, first
    InitNode* node = new InitNode;
    InitBranch* rightBranch = new InitBranch;

    // Wire rightBranch
    rightBranch->child = node;
    rightBranch->parent = branch->parent;
    rightBranch->alleles = branch->alleles; // Initialize to same as branch; will adjust
    rightBranch->nMuts = branch->nMuts;

    // Update childbranches of parent of branch
    std::erase(rightBranch->parent->childBranches, branch);
    rightBranch->parent->childBranches.push_back(rightBranch);

    // Update allBranches
    allBranches.push_back(rightBranch);

    // Rewire branch
    branch->parent = node;

    // Wire node
    node->parentBranch = rightBranch;
    node->childBranches = {branch};
    node->name = ""; // internal node

    for(size_t pos = 0; pos < branch->alleles.size(); pos++) {
        // If no mutation, nothing to do here!
        if(!isMutation(branch->alleles[pos])) {
            continue;
        }

        if(isLeft[pos]) {

          
            // If the mutation is supposed to happen on the left branch, modify the right branch
            rightBranch->alleles[pos] = getParent(branch->alleles[pos]);
            rightBranch->nMuts--;
        }else{

            // Vice versa
            branch->alleles[pos] = getChild(branch->alleles[pos]);
            branch->nMuts--;
        }
    }

    return node;

}

void attach(
    InitNode* node,
    const std::vector<uint8_t>& leafGenotype,
    InitBranch* branch,
    std::vector<InitBranch*>& allBranches
) {

    size_t L = leafGenotype.size();


    // On the existing branch, sort the mutations (2-letter words) into two categories:
    // "Left": the mutation is on the branch's child, not leafGenotype
    // "Right": the mutation is on both the branch's child and leafGenotype
    std::vector<bool> isLeft(L, false);
    int nLeft = 0;
    int nMutsOnBranch = 0; // Number of mutations on "branch"
    int nMutsOnNewBranch = 0; // Number of mutations on "newBranch"

    // Then, also record which mutations must occur to get to leafGenotype
    std::vector<uint8_t> newBranchAlleles(L);

    for(size_t pos = 0; pos < L; pos++) {
        
        if(isMutation(branch->alleles[pos])) {
            nMutsOnBranch++;
        }

        // Casework as per above comments
        // If leaf genotype is N or -, going to default to child genotype in case of mutation
        // Hence the mutation must occur on the right
        // (Can see if other way is better...)
        if(isMissing(leafGenotype[pos])) {

            if(isMutation(branch->alleles[pos])) {
                newBranchAlleles[pos] = getChild(branch->alleles[pos]);
            }else{
                newBranchAlleles[pos] = branch->alleles[pos];
            }
            
            continue;
        }

        // First, check if attachment branch has mutation or no
        if(!isMutation(branch->alleles[pos])) {

          
            // Mutation?
            if(leafGenotype[pos] != branch->alleles[pos]) {
                newBranchAlleles[pos] = makeMutation(leafGenotype[pos], branch->alleles[pos]);
                nMutsOnNewBranch++;
            }else {
                // No mutation needed - leaf matches branch
                newBranchAlleles[pos] = leafGenotype[pos];
            }
        }else{


            if(matchesChild(leafGenotype[pos], branch->alleles[pos])) {

                // Do nothing; already child allele and mutation is on the right
                newBranchAlleles[pos] = leafGenotype[pos];
                

            }else if(matchesParent(leafGenotype[pos], branch->alleles[pos])) {
                // Note that we must attach to the right of the mutation; hence, the mutation is on the left
                isLeft[pos] = true;
                nLeft++;
                newBranchAlleles[pos] = leafGenotype[pos];

            }else{
                // Very rare triallelic case: attach to left of mutation by default; hence mutation is on the right
              
                newBranchAlleles[pos] = makeMutation(leafGenotype[pos], getChild(branch->alleles[pos]));
                nMutsOnNewBranch++;

            }
        }
    }

    if(nMutsOnBranch != branch->nMuts) {
        throw std::runtime_error("Mutation count mismatch");
    }

    // Create new connecting branch
    InitBranch* newBranch = new InitBranch;
    newBranch->child = node;
    newBranch->alleles = newBranchAlleles;
    newBranch->nMuts = nMutsOnNewBranch;

    // Update allBranches
    allBranches.push_back(newBranch);

    // Wire node to newBranch
    node->parentBranch = newBranch;

    // If the following hold:
    // (1) all attachment points are left
    // (2) left node isn't a leaf
    // then attach directly to left node
    if((nLeft == 0) && (branch->child->name == "")) {
        newBranch->parent = branch->child; // Left endpoint of existing branch
        newBranch->parent->childBranches.push_back(newBranch);
        return;
    }

    // Same for right
    if(nLeft == nMutsOnBranch) {
        if(branch->parent->name != "") {
            std::cout << branch->parent->name << std::endl;
            throw std::runtime_error("Node with children shouldn't be named");
        }

        newBranch->parent = branch->parent;
        newBranch->parent->childBranches.push_back(newBranch);
        return;
    }

    InitNode* newNode = splitBranch(branch, isLeft, allBranches);

    // Wire up newBranch
    newBranch->parent = newNode;
    newBranch->parent->childBranches.push_back(newBranch);

}




void deleteInitTree(InitNode* node) {
    for(InitBranch* b : node->childBranches) {
        deleteInitTree(b->child);
        delete b;
    }
    delete node;
}

void processByLine(
    const std::string& fasta_filepath,
    const std::string& ref_filepath,
    int start_pos,
    int end_pos,
    int maxBranches,
    bool useGPU
) {
    // Genome length
    size_t L = end_pos - start_pos;

    std::string line;
    std::string refName;
    std::vector<uint8_t> refGenotype;

    // Read reference sequence
    std::ifstream ref_fasta_file(ref_filepath);
    if (!ref_fasta_file.is_open()) {
        std::cerr << "Failed to open ref_fasta_file\n";
        return;
    }

    while (std::getline(ref_fasta_file, line)) {
        if(line.empty()) continue;

        if(line[0] == '>') {
            refName = line.substr(1, line.size() - 1);
            continue;
        }

        if(start_pos < 0 || (size_t)start_pos > line.size()) {
            throw std::runtime_error("Start pos out of bounds");
        }
        if((size_t)end_pos > line.size()) {
            throw std::runtime_error("end pos too big");
        }
        if(end_pos <= start_pos) {
            throw std::runtime_error("End pos <= start_pos");
        }

        std::string letters = line.substr(start_pos, L);
        checkAlphabet(letters, {'A', 'C', 'G', 'T'});
        refGenotype = encodeString(letters);
    }
    ref_fasta_file.close();

    if(refName == "") {
        throw std::runtime_error("Invalid refName");
    }

    // PERSISTENT GPU MEMORY - allocate once for entire run
    uint8_t* d_branchAllelesFlat = nullptr;
    uint8_t* d_leafGenotypesFlat = nullptr;
    int* d_distances = nullptr;
    
    if(useGPU) {
        // Allocate for branch alleles
        cudaMallocManaged(&d_branchAllelesFlat, maxBranches * L * sizeof(uint8_t));
        cudaError_t err1 = cudaGetLastError();
        if(err1 != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA malloc failed for branches: ") + cudaGetErrorString(err1));
        }
        
        // Allocate for leaf genotypes
        cudaMallocManaged(&d_leafGenotypesFlat, maxBranches * L * sizeof(uint8_t));
        cudaError_t err2 = cudaGetLastError();
        if(err2 != cudaSuccess) {
            cudaFree(d_branchAllelesFlat);
            throw std::runtime_error(std::string("CUDA malloc failed for leaves: ") + cudaGetErrorString(err2));
        }
        
        // Allocate for distances
        cudaMallocManaged(&d_distances, maxBranches * sizeof(int));
        cudaError_t err3 = cudaGetLastError();
        if(err3 != cudaSuccess) {
            cudaFree(d_branchAllelesFlat);
            cudaFree(d_leafGenotypesFlat);
            throw std::runtime_error(std::string("CUDA malloc failed for distances: ") + cudaGetErrorString(err3));
        }
    }

    // Open main FASTA file
    std::ifstream fasta_file(fasta_filepath);
    if (!fasta_file.is_open()) {
        std::cerr << "Failed to open fasta_file\n";
        return;
    }

    // Create root node
    InitNode* rootNode = new InitNode;
    rootNode->name = "";
    rootNode->parentBranch = nullptr;
    rootNode->childBranches = {};

    // Valid tokens
    std::unordered_set<char> alphabet = {'A', 'C', 'G', 'T', 'N', '-'};

    // Track all branches
    std::vector<InitBranch*> allBranches;

    line = "";
    std::string name;

    long long total_search_time = 0;
    long long total_attach_time = 0;
    
    size_t sequenceCount = 0;

    while (std::getline(fasta_file, line)) {
        if(line.empty()) continue;

        if(line[0] == '>') {
            name = line.substr(1, line.size() - 1);
            continue;
        }

        // Check if we've hit max branches
        if(allBranches.size() >= (size_t)maxBranches) {
            std::cout << "Reached maximum number of branches (" << maxBranches << "). Stopping early.\n";
            break;
        }

        // Parse sequence
        std::string letters = line.substr(start_pos, L);
        checkAlphabet(letters, alphabet);
        std::vector<uint8_t> leafGenotype = encodeString(letters);

        // Store leaf genotype in persistent GPU memory immediately
        if(useGPU) {
            uint8_t* dest = d_leafGenotypesFlat + (sequenceCount * L);
            for(size_t j = 0; j < L; j++) {
                dest[j] = leafGenotype[j];
            }
        }

        InitNode* node = new InitNode;
        if(name == "") {
            throw std::runtime_error("Invalid sequence name");
        }
        node->name = name;
        node->childBranches = {};
        node->seq = letters;

        // FIRST BRANCH - special case
        if(allBranches.size() == 0) {
            InitBranch* branch = new InitBranch;
            branch->parent = rootNode;
            branch->child = node;
            branch->alleles = std::vector<uint8_t>(L, 0);
            branch->nMuts = 0;

            for(size_t pos = 0; pos < L; pos++) {
                if(isMissing(leafGenotype[pos])) {
                    branch->alleles[pos] = refGenotype[pos];
                    continue;
                }
                if(leafGenotype[pos] != refGenotype[pos]) {
                    branch->alleles[pos] = makeMutation(leafGenotype[pos], refGenotype[pos]);
                    branch->nMuts++;
                } else {
                    branch->alleles[pos] = leafGenotype[pos];
                }
            }

            node->parentBranch = branch;
            rootNode->childBranches.push_back(branch);
            allBranches.push_back(branch);

            // Copy first branch to GPU
            if(useGPU) {
                for(size_t j = 0; j < L; j++) {
                    d_branchAllelesFlat[j] = branch->alleles[j];
                }
            }

            sequenceCount++;
            continue;
        }

        // SEARCH FOR BEST ATTACHMENT BRANCH
        auto search_start = std::chrono::high_resolution_clock::now();

        size_t numBranches = allBranches.size();
        int bestDist = std::numeric_limits<int>::max();
        InitBranch* bestBranch = nullptr;

        if(useGPU) {
            // ============ GPU PATH ============
            
            // Get pointer to current leaf genotype (already in GPU memory!)
            uint8_t* d_leafGenotype = d_leafGenotypesFlat + (sequenceCount * L);

            // Call GPU kernel with persistent memory - no allocations!
            computeDistancesGPU(d_leafGenotype, d_branchAllelesFlat, d_distances, numBranches, L);

            // Find minimum distance on CPU
            for(size_t i = 0; i < numBranches; i++) {
                if(d_distances[i] < bestDist) {
                    bestDist = d_distances[i];
                    bestBranch = allBranches[i];
                }
                if(d_distances[i] == 0) break;
            }

        } else {
            // ============ CPU PATH ============
            for(InitBranch* b : allBranches) {
                int dist = evoDistance(leafGenotype, b->alleles);
                if(dist < bestDist) {
                    bestDist = dist;
                    bestBranch = b;
                }
                if(dist == 0) break;
            }
        }

        auto search_end = std::chrono::high_resolution_clock::now();
        total_search_time += std::chrono::duration_cast<std::chrono::microseconds>(
            search_end - search_start).count();

        // Find index of branch being attached to
        size_t attachBranchIdx = SIZE_MAX;
        for(size_t i = 0; i < allBranches.size(); i++) {
            if(allBranches[i] == bestBranch) {
                attachBranchIdx = i;
                break;
            }
        }

        size_t numBranchesBefore = allBranches.size();

        // Attach to tree
        auto attach_start = std::chrono::high_resolution_clock::now();
        attach(node, leafGenotype, bestBranch, allBranches);
        auto attach_end = std::chrono::high_resolution_clock::now();
        total_attach_time += std::chrono::duration_cast<std::chrono::microseconds>(
            attach_end - attach_start).count();

        // Update GPU memory for changed branches
        if(useGPU) {
            // Check if branch was split (creates 2 new branches instead of 1)
            bool branchWasSplit = (allBranches.size() > numBranchesBefore + 1);
            
            if(branchWasSplit && attachBranchIdx != SIZE_MAX) {
                // Original branch was modified - update it on GPU
                uint8_t* dest = d_branchAllelesFlat + (attachBranchIdx * L);
                for(size_t j = 0; j < L; j++) {
                    dest[j] = allBranches[attachBranchIdx]->alleles[j];
                }
            }
            
            // Copy new branches to GPU (1-2 new branches)
            for(size_t i = numBranchesBefore; i < allBranches.size(); i++) {
                uint8_t* dest = d_branchAllelesFlat + (i * L);
                for(size_t j = 0; j < L; j++) {
                    dest[j] = allBranches[i]->alleles[j];
                }
            }
        }
        
        sequenceCount++;
    }
    
    fasta_file.close();

    // Print performance stats
    std::cout << "\n=== Performance Stats ===\n";
    std::cout << "Mode: " << (useGPU ? "GPU" : "CPU") << "\n";
    std::cout << "Sequences processed: " << sequenceCount << "\n";
    std::cout << "Total search time: " << total_search_time / 1000.0 << " ms\n";
    std::cout << "Total attach time: " << total_attach_time / 1000.0 << " ms\n";

    int totMutations = 0;
    for(InitBranch* b : allBranches) {
        totMutations += b->nMuts;
    }
    std::cout << "Total Mutations: " << totMutations << std::endl;

    // Free GPU memory once at the end
    if(useGPU) {
        cudaFree(d_branchAllelesFlat);
        cudaFree(d_leafGenotypesFlat);
        cudaFree(d_distances);
    }

    // Validate tree structure
    checkTree(refGenotype, rootNode);

    // Cleanup tree
    deleteInitTree(rootNode);
}

int main(int argc, char* argv[]) {
    if(argc != 7) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <ref_file> <start_pos> <end_pos> <max_branches> <use_gpu>\n";
        std::cerr << "  use_gpu: 0 for CPU, 1 for GPU\n";
        std::cerr << "Example: " << argv[0] << " data.fasta ref.fasta 0 30000 10000 1\n";
        return 1;
    }

    std::string fasta_filepath = argv[1];
    std::string ref_filepath = argv[2];
    int start_pos = std::stoi(argv[3]);
    int end_pos = std::stoi(argv[4]);
    int maxBranches = std::stoi(argv[5]);
    bool useGPU = (std::stoi(argv[6]) != 0);

    std::cout << "Running in " << (useGPU ? "GPU" : "CPU") << " mode\n";

    processByLine(fasta_filepath, ref_filepath, start_pos, end_pos, maxBranches, useGPU);

    std::cout << "Done!" << std::endl;

    return 0;
}