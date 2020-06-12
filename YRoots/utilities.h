//
//  utilities.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef utilities_h
#define utilities_h

#include <vector>
#include <string>
#include <math.h>

#define _USE_MATH_DEFINES
#include <cmath>

//TODO: Make unit tests for these functions!

std::vector<std::string> split(const std::string& string, const std::string& delimiter) {
    std::vector<std::string> results;
    size_t last = 0;
    size_t next = 0;
    while ((next = string.find(delimiter, last)) != std::string::npos) {
        results.push_back(string.substr(last, next-last));
        last = next + 1;
    }
    results.push_back(string.substr(last));
    return results;
}

bool is_number(const std::string& s) {
    return( strspn( s.c_str(), "-.0123456789" ) == s.size() );
}

//TODO: Check for overflow!
size_t power(size_t base, size_t exponent) {
    size_t result = 1;
    for(size_t i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}

enum Dimension {
    NDim,
    One,
    Two,
    Three,
    Four
};

struct Interval {
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    
    Interval() {}
};

#endif /* utilities_h */
