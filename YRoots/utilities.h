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
#include <iostream>

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

struct SubdivisionParameters {
    double relApproxTol = 1.e-15;
    double absApproxTol = 1.e-12;
    double maxConditionNumber = 1e5;
    double goodZerosFactor = 100;
    double minGoodZerosTol = 1e-5;
    bool checkEvaluationError = true;
    size_t checkEvaluationFrequency = 1;
    size_t approximationDegree = 3;
    size_t targetDegree = 10;
    size_t maxLevel = 999;
    bool returnPotentials = false;
    std::string method = "svd";
    double targetTol = 1.e-15;
    bool useTargetTol = true;
};

struct SolveParameters {
    Interval interval;
    size_t currentLevel;
    
    SolveParameters() : currentLevel(0) {
        
    }
    
    SolveParameters(Interval _interval, size_t _currentLevel) :
    interval(_interval),
    currentLevel(_currentLevel)
    {
        
    }
};

#endif /* utilities_h */
