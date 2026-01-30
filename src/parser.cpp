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
#include <cstring>
#include "parser.h"

double date_to_decimal_year(int y, int m, int d)
{
    static const int mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31};

    auto is_leap = [](int y) {
        return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
    };

    int doy = d - 1;
    for (int i = 0; i < m - 1; ++i)
        doy += mdays[i];

    if (m > 2 && is_leap(y))
        doy += 1;

    int days_in_year = is_leap(y) ? 366 : 365;

    return y + double(doy) / days_in_year;
}

std::unordered_map<std::string, double> parseMetadata(const std::string& metadata_filepath) {
    std::cout << "Reading entire metadata file into memory..." << std::flush;
    std::ifstream metadata_file(metadata_filepath);
    std::stringstream buffer;
    buffer << metadata_file.rdbuf();
    std::string metadata_content = buffer.str();
    metadata_file.close();
    std::cout << " OK (" << metadata_content.size() / 1e6 << " MB)\n";

    std::unordered_map<std::string, double> out;
    out.reserve(100000);  // Pre-allocate to avoid rehashing

    const char* data = metadata_content.c_str();
    const char* end = data + metadata_content.size();
    const char* ptr = data;

    size_t lines_processed = 0;
    size_t entries_parsed = 0;

    std::cout << "Parsing metadata entries...\n";

    while (ptr < end) {
        // Find start of line
        const char* line_start = ptr;
        
        // Find end of line
        const char* line_end = (const char*)memchr(ptr, '\n', end - ptr);
        if (!line_end) line_end = end;
        
        const char* accession_start = nullptr;
        const char* accession_end = nullptr;
        const char* date_start = nullptr;
        const char* date_end = nullptr;

        // Search for "accession":"
        const char* search_ptr = line_start;
        while (search_ptr < line_end - 13) {
            if (memcmp(search_ptr, "\"accession\":\"", 13) == 0) {
                accession_start = search_ptr + 13;
                // Find closing quote
                accession_end = (const char*)memchr(accession_start, '\"', line_end - accession_start);
                break;
            }
            search_ptr++;
        }

        // Search for "collectionDate":"
        search_ptr = line_start;
        while (search_ptr < line_end - 18) {
            if (memcmp(search_ptr, "\"collectionDate\":\"", 18) == 0) {
                date_start = search_ptr + 18;
                // Find closing quote
                date_end = (const char*)memchr(date_start, '\"', line_end - date_start);
                break;
            }
            search_ptr++;
        }

        // If we found both fields, parse the date
        if (accession_start && accession_end && date_start && date_end) {
            // Parse date in format YYYY-MM-DD
            // Count dashes to ensure it's complete
            int dash_count = 0;
            for (const char* p = date_start; p < date_end; p++) {
                if (*p == '-') dash_count++;
            }

            if (dash_count == 2) {
                // Fast integer parsing
                int y = 0, m = 0, d = 0;
                const char* p = date_start;
                
                // Parse year (4 digits)
                while (p < date_end && *p != '-') {
                    y = y * 10 + (*p - '0');
                    p++;
                }
                p++; // skip dash
                
                // Parse month (1-2 digits)
                while (p < date_end && *p != '-') {
                    m = m * 10 + (*p - '0');
                    p++;
                }
                p++; // skip dash
                
                // Parse day (1-2 digits)
                while (p < date_end) {
                    d = d * 10 + (*p - '0');
                    p++;
                }

                double decimal_date = date_to_decimal_year(y, m, d);
                
                // Create string from accession
                std::string seq_name(accession_start, accession_end - accession_start);
                out[seq_name] = decimal_date;
                
                entries_parsed++;
                
                // Log every 10000 entries
                if (entries_parsed % 10000 == 0) {
                    std::cout << "  Parsed " << entries_parsed << " entries...\n" << std::flush;
                }
            }
        }

        lines_processed++;
        
        // Move to next line
        ptr = line_end + 1;
    }

    std::cout << "Parsed " << out.size() << " metadata entries from " 
              << lines_processed << " lines\n";
    return out;
}