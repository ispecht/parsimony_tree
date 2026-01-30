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

    for(size_t pos = 0; pos < seq.size(); pos++) {
        out.push_back(encodeBase(seq[pos]));
    }

    return out;
}


double dPois(int x, double lambda, bool log) {
    // Input validation
    if (x < 0) {
        throw std::invalid_argument("x must be non-negative");
    }

    if(lambda == 0 && x == 0){
        if(log){
            return 0.0;
        }else{
            return 1.0;
        }
    }



    if (lambda <= 0) {
        throw std::invalid_argument("lambda must be positive");
    }
    
    // For x = 0, the formula simplifies to e^(-λ)
    if (x == 0) {
        if(log) {
            return (-lambda);
        }else{
            return std::exp(-lambda);
        }
    }
    
    // Use logarithms to avoid overflow for large values
    // log(P(X = x)) = x * log(λ) - λ - log(x!)
    // We compute log(x!) using lgamma(x + 1)
    double log_pmf = x * std::log(lambda) - lambda - std::lgamma(x + 1);

    if(log) {
        return log_pmf;
    }else{
        return std::exp(log_pmf);
    }
}



double rUnif(double a, double b, std::mt19937& rng) {
    // Input validation
    if (a > b) {
        throw std::invalid_argument("Parameter 'a' must be less than 'b'");
    }

    // Set up random number generation
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}