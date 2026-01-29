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
#include "helpers.h"

uint8_t encodeBase(char c) {
    if(c == 'A') {
        return 1;
    }else if(c == 'C') {
        return 2;
    }else if(c == 'G') {
        return 4;
    }else if (c == 'T') {
        return 8;
    }else{
        return 0;
    }
}

char decodeBase(uint8_t encoding) {
    if(encoding == 1) {
        return 'A';
    }else if(encoding == 2) {
        return 'C';
    }else if(encoding == 4) {
        return 'G';
    }else if(encoding == 8) {
        return 'T';
    }else{
        return 'N';
    }
}

// First is child allele
// Second is parent allele
uint8_t makeMutation(uint8_t first, uint8_t second) {
    return (second << 4) | first;
}

bool matchesChild(uint8_t allele, uint8_t mutation) {
    return (mutation & 0x0F) == allele;
}

// Check if encoded allele matches the parent/second allele (high nibble)
bool matchesParent(uint8_t allele, uint8_t mutation) {
    return (mutation >> 4) == allele;
}

uint8_t getChild(uint8_t mutation) {
    return (mutation & 0x0F);
}

// Check if encoded allele matches the parent/second allele (high nibble)
uint8_t getParent(uint8_t mutation) {
    return (mutation >> 4);
}

bool isMutation(uint8_t allele) {
    return (allele >> 4);
}

bool isMissing(uint8_t allele) {
    return (allele == 0);
}

std::vector<uint8_t> encodeString(std::string seq) {
    std::vector<uint8_t> out;

    for(int pos = 0; pos < seq.size(); pos++) {
        out.push_back(encodeBase(seq[pos]));
    }

    return out;
}