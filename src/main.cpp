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
#include "parser.h"

// Add this structure to hold sequence metadata
struct SequenceMetadata {
    std::string name;
    double date;
    std::streampos file_position;  // Position in file where sequence starts
    size_t sequence_length;
};


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
    const std::string& metadata_filepath,
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

    // Read entire file into memory first

    // What we actually need to do:
    // - Get every sequence in the fasta (vector)
    // - Get the name of each sequence
    // - Based on metadata, get time of each sequence

    // Parse metadata
    std::unordered_map<std::string, double> name_to_date = parseMetadata(metadata_filepath);

    // PHASE 1: Collect sequence metadata without loading sequences
    std::cout << "Phase 1: Scanning FASTA for sequence metadata..." << std::flush;
    std::vector<SequenceMetadata> sequence_metadata;

    std::ifstream fasta_scan(fasta_filepath);
    if (!fasta_scan.is_open()) {
        throw std::runtime_error("Failed to open FASTA file");
    }

    std::string current_name;  // Only declare new variables here
    std::streampos name_position = 0;

    while (std::getline(fasta_scan, line)) {  // Use existing 'line' variable
        if(line.empty()) continue;
        
        if(line[0] == '>') {
            std::string full_name = line.substr(1);
            size_t space_pos = full_name.find(' ');
            current_name = (space_pos != std::string::npos) ? full_name.substr(0, space_pos) : full_name;
            name_position = fasta_scan.tellg();
        } else if (!current_name.empty()) {
            auto it = name_to_date.find(current_name);
            if(it != name_to_date.end()) {
                SequenceMetadata meta;
                meta.name = current_name;
                meta.date = it->second;
                meta.file_position = name_position;
                meta.sequence_length = line.length();
                sequence_metadata.push_back(meta);
            }
            current_name.clear();
        }
    }
    fasta_scan.close();

    std::cout << " Found " << sequence_metadata.size() << " sequences\n";

    // PHASE 2: Sort by date
    std::cout << "Phase 2: Sorting sequences by date..." << std::flush;
    std::sort(sequence_metadata.begin(), sequence_metadata.end(),
            [](const SequenceMetadata& a, const SequenceMetadata& b) {
                return a.date < b.date;
            });
    std::cout << " OK\n";

    // PERSISTENT GPU MEMORY - allocate once for entire run
    uint8_t* d_branchAllelesFlat = nullptr;
    uint8_t* d_leafGenotype = nullptr;
    int* d_distances = nullptr;
    
    if(useGPU) {

        // Check if GPU is available
        int deviceCount = 0;
        cudaError_t err = cudaGetDeviceCount(&deviceCount);
        if(err != cudaSuccess || deviceCount == 0) {
            throw std::runtime_error("No CUDA-capable GPU found!");
        }
        
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        std::cout << "Using GPU: " << prop.name << "\n";
        std::cout << "GPU memory: " << prop.totalGlobalMem / 1e9 << " GB\n";

        size_t branchBytes = maxBranches * L * sizeof(uint8_t);
        size_t distBytes = maxBranches * sizeof(int);
        size_t totalBytes = branchBytes + distBytes;  // Remove leafBytes
        
        std::cout << "Allocating GPU memory:\n";
        std::cout << "  Branches: " << branchBytes / 1e6 << " MB\n";
        std::cout << "  Distances: " << distBytes / 1e6 << " MB\n";
        std::cout << "  Total: " << totalBytes / 1e6 << " MB\n";
        
        // Allocate for branch alleles
        std::cout << "Allocating d_branchAllelesFlat..." << std::flush;
        cudaError_t err1 = cudaMallocManaged(&d_branchAllelesFlat, branchBytes);
        if(err1 != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA malloc failed for branches: ") + cudaGetErrorString(err1));
        }
        std::cout << " OK\n";
        
        // Allocate for distances
        std::cout << "Allocating d_distances..." << std::flush;
        cudaError_t err3 = cudaMallocManaged(&d_distances, distBytes);
        if(err3 != cudaSuccess) {
            cudaFree(d_branchAllelesFlat);
            throw std::runtime_error(std::string("CUDA malloc failed for distances: ") + cudaGetErrorString(err3));
        }
        std::cout << " OK\n";

        // Allocate single reusable leaf buffer
        std::cout << "Allocating d_leafGenotype (single)..." << std::flush;
        cudaError_t err_leaf = cudaMallocManaged(&d_leafGenotype, L * sizeof(uint8_t));
        if(err_leaf != cudaSuccess) {
            cudaFree(d_branchAllelesFlat);
            cudaFree(d_distances);
            throw std::runtime_error(std::string("CUDA malloc failed for leaf: ") + cudaGetErrorString(err_leaf));
        }
        std::cout << " OK\n";

    }


    // Create root node
    InitNode* rootNode = new InitNode;
    rootNode->name = "";
    rootNode->parentBranch = nullptr;
    rootNode->childBranches = {};

    // Valid tokens
    std::unordered_set<char> alphabet = {'A', 'C', 'G', 'T'};

    // Track all branches
    std::vector<InitBranch*> allBranches;

    line = "";
    std::string name;

    long long total_search_time = 0;
    long long total_attach_time = 0;
    long long total_gpu_update_time = 0;  // NEW

    
    size_t sequenceCount = 0;

    
    
    // PHASE 3: Process sequences in chronological order
    std::cout << "Phase 3: Processing sequences in chronological order...\n";

    // Reopen file for reading sequences
    std::ifstream fasta_file(fasta_filepath);
    if (!fasta_file.is_open()) {
        throw std::runtime_error("Failed to reopen FASTA file for processing");
    }
    std::string sequence_line;

    for(const auto& meta : sequence_metadata) {
        if(sequenceCount % 1000 == 0) {
            std::cout << "Processing sequence " << sequenceCount 
                    << " (date: " << meta.date << ")\n";
        }
        
        // Check branch limit
        if(allBranches.size() + 2 >= (size_t)maxBranches) {
            std::cout << "Reached maximum number of branches (" << maxBranches 
                    << "). Stopping early.\n";
            break;
        }
        
        // Seek to sequence position and read
        fasta_file.seekg(meta.file_position);
        if(!std::getline(fasta_file, sequence_line)) {
            std::cerr << "Warning: Failed to read sequence for " << meta.name << "\n";
            continue;
        }
        
        // Parse sequence
        std::string letters = sequence_line.substr(start_pos, L);
        std::vector<uint8_t> leafGenotype = encodeString(letters);
        
        // Copy to GPU
        if(useGPU) {
            cudaMemcpy(d_leafGenotype, leafGenotype.data(), 
                      L * sizeof(uint8_t), cudaMemcpyHostToDevice);
        }
        
        InitNode* node = new InitNode;
        node->name = meta.name;
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
                //std::cout << "  Copying first branch to GPU..." << std::flush;
                for(size_t j = 0; j < L; j++) {
                    d_branchAllelesFlat[j] = branch->alleles[j];
                }
                //std::cout << " OK\n" << std::flush;
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
            
            //std::cout << "  Calling GPU kernel (numBranches=" << numBranches << ")..." << std::flush;
            computeDistancesGPU(d_leafGenotype, d_branchAllelesFlat, d_distances, numBranches, L);
            //std::cout << " OK\n" << std::flush;

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

        if(useGPU) {
            auto gpu_update_start = std::chrono::high_resolution_clock::now();
            
            bool branchWasSplit = (allBranches.size() > numBranchesBefore + 1);
            
            if(branchWasSplit && attachBranchIdx != SIZE_MAX) {
                uint8_t* dest = d_branchAllelesFlat + (attachBranchIdx * L);
                cudaMemcpy(dest, allBranches[attachBranchIdx]->alleles.data(), 
                        L * sizeof(uint8_t), cudaMemcpyHostToDevice);
            }
            
            for(size_t i = numBranchesBefore; i < allBranches.size(); i++) {
                uint8_t* dest = d_branchAllelesFlat + (i * L);
                cudaMemcpy(dest, allBranches[i]->alleles.data(), 
                        L * sizeof(uint8_t), cudaMemcpyHostToDevice);
            }
            
            auto gpu_update_end = std::chrono::high_resolution_clock::now();
            total_gpu_update_time += std::chrono::duration_cast<std::chrono::microseconds>(
                gpu_update_end - gpu_update_start).count();
        }
        
        sequenceCount++;
    }
    
    fasta_file.close();

    std::cout << "\n=== Performance Stats ===\n";
    std::cout << "Mode: " << (useGPU ? "GPU" : "CPU") << "\n";
    std::cout << "Sequences processed: " << sequenceCount << "\n";
    std::cout << "Total search time: " << total_search_time / 1000.0 << " ms\n";
    std::cout << "Total attach time: " << total_attach_time / 1000.0 << " ms\n";
    if(useGPU) {
        std::cout << "Total GPU update time: " << total_gpu_update_time / 1000.0 << " ms\n";
    }

    int totMutations = 0;
    for(InitBranch* b : allBranches) {
        totMutations += b->nMuts;
    }
    std::cout << "Total Mutations: " << totMutations << std::endl;

    // Free GPU memory once at the end
    if(useGPU) {
        cudaFree(d_branchAllelesFlat);
        cudaFree(d_leafGenotype);  // CHANGED: single buffer
        cudaFree(d_distances);
    }

    bool SAFETY_MODE = false;
    // Validate tree structure
    if(SAFETY_MODE) {
        checkTree(refGenotype, rootNode);

        std::cout << "Mutation counts: " << std::endl;

        for(InitBranch* b : allBranches) {
            std::cout << b->nMuts << std::endl;
        }
    }

    // Cleanup tree
    deleteInitTree(rootNode);
}

int main(int argc, char* argv[]) {
    if(argc != 8) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <ref_file> <metadata_file> <start_pos> <end_pos> <max_branches> <use_gpu>\n";
        std::cerr << "  use_gpu: 0 for CPU, 1 for GPU\n";
        std::cerr << "Example: " << argv[0] << " data.fasta ref.fasta 0 30000 10000 1\n";
        return 1;
    }

    std::string fasta_filepath = argv[1];
    std::string ref_filepath = argv[2];
    std::string metadata_filepath = argv[3];
    int start_pos = std::stoi(argv[4]);
    int end_pos = std::stoi(argv[5]);
    int maxBranches = std::stoi(argv[6]);
    bool useGPU = (std::stoi(argv[7]) != 0);

    std::cout << "Running in " << (useGPU ? "GPU" : "CPU") << " mode\n";

    try {
        processByLine(fasta_filepath, ref_filepath, metadata_filepath, start_pos, end_pos, maxBranches, useGPU);
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Done!" << std::endl;

    return 0;
}