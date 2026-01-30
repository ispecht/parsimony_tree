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

std::unordered_map<std::string, double> parseMetadata(const std::string& metadata_filepath);