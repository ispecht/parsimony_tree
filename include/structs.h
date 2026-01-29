#pragma once

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
    std::vector<uint8_t> alleles; // Allowed alleles at each position
    int nMuts;
};