//
//  utilities.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef utilities_h
#define utilities_h

#include <algorithm>
#include <cctype>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <cmath>
#include <stddef.h>
#include "../Libraries/Eigen/Core"
#include "../Utilities/macros.hpp"

//TODO: Make unit tests for these functions!

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

std::string formatTimePretty(double nanoseconds) {
    static constexpr double thousand = 1000;
    static constexpr double million = thousand*thousand;
    static constexpr double billion = million*thousand;

    if(nanoseconds < thousand) {
        return std::to_string(nanoseconds) + "ns";
    }
    else if(nanoseconds < million) {
        return std::to_string(nanoseconds/thousand) + "us";
    }
    else if(nanoseconds < billion) {
        return std::to_string(nanoseconds/million) + "ms";
    }
    else {
        return std::to_string(nanoseconds/billion) + "s";
    }
}

void printMatrix(const Eigen::MatrixXd& matrix) {
    for(size_t i = 0; i < static_cast<size_t>(matrix.rows()); i++) {
        for(size_t j = 0; j < static_cast<size_t>(matrix.cols()); j++) {
            std::cout<<matrix(i,j)<<"\t";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}

void printVector(const Eigen::VectorXd& vector) {
    for(size_t i = 0; i < static_cast<size_t>(vector.rows()); i++) {
        std::cout<<vector(i)<<"\t";
    }
    std::cout<<"\n\n";
}

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

void replaceStringInPlace(std::string& subject, const std::string& search, const std::string& replace) {
    //Copied from https://stackoverflow.com/questions/5878775/how-to-find-and-replace-string/5878802
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}

bool is_number(const std::string& s) {
    return( strspn( s.c_str(), "-.0123456789" ) == s.size() );
}

bool isDigit(const char c) {
    //Return if c is in "-.0123456789"
    return c >= 45 && c <= 57 && c != 47;
}

bool isNumericDigit(const char c) {
    //Return if c is in "0123456789"
    return c >= 48 && c <= 57;
}


inline void printWarning(const std::string& warningMessage) {
    //Add a newline before and after to not interfere with the percent tracking.
    std::cout<<"\n"<<warningMessage<<"\n";
}

inline void printAndThrowRuntimeError(const std::string& errorMessage) {
    //std::cout<<errorMessage<<"\n"; //I don't really need to print it, the runtime error prints it.
    throw std::runtime_error(errorMessage);
}

std::string toLowerSubstring(std::string _input, size_t _start, size_t _numChars) {
    std::string output;
    if(unlikely(_start + _numChars >= _input.length())) {
        printAndThrowRuntimeError("Invalid Subtring!");
    }
    else {
        for(size_t i = _start; i < _start + _numChars; i++) {
            output += std::tolower(_input[i]);
        }
    }
    return output;
}

template<typename T1, typename T2>
inline T1 power(T1 base, T2 exponent) {
    return std::pow(base, exponent);
}

template<typename T>
T power(T base, size_t exponent) {
    //Base Cases
    if (exponent == 0) {
         return 1;
    }
    if (exponent == 1) {
        return base;
    }

    //Log(n) complexity recursion
    T temp = power(base, exponent/2);
    if ((exponent&1) == 0) { //Fast Check if exponent is even
        return temp * temp;
    }
    else {
        return base * temp * temp;
    }
}

template<typename T>
T chebPower(T base, size_t exponent) {
    //Use the identities T_2n(x) = 2Tn(x)^2 - T_0(x)
    //                   T_2n+1(x) = 2Tn+1(x)Tn(x) - T_1(x)
    //These will recompute most calculations after the first use, so do bottom up
    //to get the beginning vals and then use it once at the end
    
    //Base Cases
    if (exponent == 0) {
         return 1;
    }
    if (exponent == 1) {
        return base;
    }
    
    //Variables for calculating linear part
    double twoBase = 2*base;
    double temp;
    double temp0 = 1;
    double temp1 = base;

    //Log(n) complexity recursion
    if ((exponent&1) == 0) { //Fast Check if exponent is even
        size_t toCompute = exponent / 2 - 1;
        //Temp1 is now T_(toCompute + 1)
        for(size_t i = 0; i < toCompute; i++) {
            temp = temp1;
            temp1 = twoBase * temp1 - temp0;
            temp0 = temp;
        }
        //Use the identity T_2n(x) = 2Tn(x)^2 - T_0(x)
        return 2 * temp1 * temp1 - 1;
    }
    else {
        size_t toCompute = exponent / 2;
        //Temp1 is now T_(toCompute + 1)
        for(size_t i = 0; i < toCompute; i++) {
            temp = temp1;
            temp1 = twoBase * temp1 - temp0;
            temp0 = temp;
        }
        //Use the identity T_2n+1(x) = 2Tn+1(x)Tn(x) - T_1(x)
        return 2 * temp1 * temp0 - base;
    }
}


enum Dimension {
    NDim,
    One,
    Two,
    Three,
    Four,
    Five
};

enum SolveMethod {
    ConstantTermCheck = 0,
    QuadraticCheck = 1,
    BoundingInterval = 2,
    LinearSolve = 3,
    SpectralSolve = 4,
    TooDeep = 5
};

struct Interval {
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    double              area;
    bool                areaFound;
    
    std::string toString() const {
        std::string result = "";
        for(size_t dim = 0; dim < lowerBounds.size(); dim++) {
            if(dim != 0) {
                result += ",";
            }
            result += "[" + std::to_string(lowerBounds[dim]) + "," + std::to_string(upperBounds[dim]) + "]";
        }
        if(lowerBounds.size() > 1) {
            result = "[" + result + "]";
        }

        return result;
    }
    
    double getArea() {
        if(areaFound) {
            return area;
        }
        
        area = 1.0;
        for(size_t i = 0; i < lowerBounds.size(); i++) {
            area *= (upperBounds[i] - lowerBounds[i]);
        }
        areaFound = true;
        return area;
    }
    
    void setArea(double _area) {
        area = _area;
        areaFound = true;
    }
    
    inline void clear() {
        areaFound = false;
    }
    
    Interval() : areaFound(false) {}
};

void projectInterval(Interval& resultInterval, const Interval& currentInterval, const Interval& projectOntoInterval) {
    for(size_t i = 0; i < currentInterval.upperBounds.size(); i++) {
        double temp1 = currentInterval.upperBounds[i] - currentInterval.lowerBounds[i];
        double temp2 = currentInterval.upperBounds[i] + currentInterval.lowerBounds[i];
        resultInterval.lowerBounds[i] = ((temp1 * projectOntoInterval.lowerBounds[i]) + temp2) /2.0;
        resultInterval.upperBounds[i] = ((temp1 * projectOntoInterval.upperBounds[i]) + temp2) /2.0;
    }
}

struct SubdivisionParameters {
    //Approx Tols
    double relApproxTol = 1.e-10; //TODO: Test to figure out the best defaults.
    double absApproxTol = 1.e-10;
    double targetTol = 1e-15;
    bool check_eval_error = true;
    //Good Zeros
    double goodZerosFactor = 100;
    double minGoodZerosTol = 1e-5;
    //Degrees
    size_t approximationDegree = 20;
    size_t targetDegree = 1;
    size_t maxLevel = 50;
};

struct GeneralParameters {
    bool trackIntervals = false;
    bool trackProgress = true;
    bool useTimer = false;
    size_t numThreads = 1;
};

struct SolveParameters {
    Interval interval;
    size_t currentLevel;
    std::vector<size_t> goodDegrees;
    
    SolveParameters() : currentLevel(0) {}
        
    inline void clear() {
        interval.clear();
    }
};


#endif /* utilities_h */
