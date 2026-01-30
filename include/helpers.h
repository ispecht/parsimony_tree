#pragma once


uint8_t encodeBase(char c);

char decodeBase(uint8_t encoding);

// First is child allele
// Second is parent allele
uint8_t makeMutation(uint8_t first, uint8_t second);

bool matchesChild(uint8_t allele, uint8_t mutation);

// Check if encoded allele matches the parent/second allele (high nibble)
bool matchesParent(uint8_t allele, uint8_t mutation);

uint8_t getChild(uint8_t mutation);

// Check if encoded allele matches the parent/second allele (high nibble)
uint8_t getParent(uint8_t mutation);

bool isMutation(uint8_t allele);

bool isMissing(uint8_t allele);

std::vector<uint8_t> encodeString(std::string seq);


double dPois(int x, double lambda, bool log);

double rUnif(double a, double b, std::mt19937& rng);