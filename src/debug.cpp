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
#include "debug.h"
#include "structs.h"
#include "helpers.h"


void checkAlphabet(const std::string& leafGenotype, const std::unordered_set<char>& alphabet) {

    for(size_t i = 0; i < leafGenotype.size(); i++) {
        if(!alphabet.contains(leafGenotype[i])) {
            throw std::runtime_error("Invalid character");
        }
    }

}

void checkTreeAtPosRec(
    uint8_t allele,
    size_t pos, 
    InitNode* node
) {

    if(node->seq != "") {
        char letter = decodeBase(allele);
        if(letter == 'N') {
            throw std::runtime_error("Decoding should never yield N");
        }
        
        if((node->seq[pos] != 'N') && (node->seq[pos] != '-')){
            if(letter != node->seq[pos]){
                throw std::runtime_error("Mismatch at leaf");
            } 
        }
    }else{

        for(InitBranch* b : node->childBranches) {

            if(isMutation(b->alleles[pos])) {

                if(getParent(b->alleles[pos]) != allele) {
                    throw std::runtime_error("Mutation mismatch");
                }
                
                // Recurse
                checkTreeAtPosRec(getChild(b->alleles[pos]), pos, b->child);
            }else{

                if(b->alleles[pos] != allele) {
                    std::cout << "Position: " << pos << std::endl;
                    std::cout << "Branch allele: " << decodeBase(b->alleles[pos]) << std::endl;
                    std::cout << "Expected allele: " << decodeBase(allele) << std::endl;
                    std::cout << "Node name: " << b->child->name << std::endl;
                    throw std::runtime_error("Mismatch between current allele and allele on branch");
                }

                checkTreeAtPosRec(allele, pos, b->child);
            }
        }
    }
}

void checkTree(
    const std::vector<uint8_t>& refGenotype,
    InitNode* rootNode
) {
    for(size_t pos = 0; pos < refGenotype.size(); pos++) {
        checkTreeAtPosRec(refGenotype[pos], pos, rootNode);
    }
}