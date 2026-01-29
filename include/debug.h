#pragma once

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
#include "structs.h"


void checkAlphabet(const std::string& leafGenotype, const std::unordered_set<char>& alphabet);

void checkTree(
    const std::vector<uint8_t>& refGenotype,
    InitNode* rootNode
);