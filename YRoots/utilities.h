//
//  utilities.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef utilities_h
#define utilities_h

#include "macros.h"

#include <algorithm>
#include <cctype>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <cmath>
#include <Eigen/Core>

//TODO: Make unit tests for these functions!

void printMatrix(Eigen::MatrixXd& matrix) {
    for(size_t i = 0; i < matrix.rows(); i++) {
        for(size_t j = 0; j < matrix.cols(); j++) {
            std::cout<<matrix(i,j)<<"\t";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}

void printVector(Eigen::VectorXd& vector) {
    for(size_t i = 0; i < vector.rows(); i++) {
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

bool is_number(const std::string& s) {
    return( strspn( s.c_str(), "-.0123456789" ) == s.size() );
}

bool isDigit(const char c) {
    //Return if c is in "-.0123456789"
    return c >= 45 && c <= 57 && c != 47;
}

bool isNumericDigit(const char c) {
    //Return if c is in "123456789"
    return c >= 49 && c <= 57;
}

void printAndThrowRuntimeError(const std::string& errorMessage) {
    std::cout<<errorMessage<<"\n";
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
    
    void clear() {
        areaFound = false;
    }
    
    Interval() {}
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
    double relApproxTol = 1.e-10;
    double absApproxTol = 1.e-10;
    double goodZerosFactor = 100;
    double minGoodZerosTol = 1e-5;
    size_t approximationDegree = 20;
    size_t targetDegree = 1;
    size_t maxLevel = 999;
    bool trackIntervals = true;
};

struct SolveParameters {
    Interval interval;
    size_t currentLevel;
    std::vector<size_t> goodDegrees;
    
    SolveParameters() : currentLevel(0) {
        
    }
    
    SolveParameters(Interval _interval, size_t _currentLevel, std::vector<size_t> _goodDegrees) :
    interval(_interval),
    currentLevel(_currentLevel),
    goodDegrees(_goodDegrees)
    {
        
    }
};





struct TimingDetails {
    //TODO: Have this hold a vector of the individual run times, so we can print the median run time, highest 10%, lowest 10% etc
    
    std::string                                                 name;
    bool                                                        indexClaimed;
    size_t                                                      runCount;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double                                                      runTimesNanos;
    
    TimingDetails() : indexClaimed(false), runCount(0), runTimesNanos(0.0) {}
    
    void claim(std::string _name) {
        if(indexClaimed && _name != name) {
            printAndThrowRuntimeError(_name + " is trying to claim a timing detail that is already claimed by " + name);
        }
        else {
            indexClaimed = true;
            name = _name;
        }
    }
    
    void clearClaim() {
        name = "";
        indexClaimed = false;
        runCount = 0;
        runTimesNanos = 0;
    }
    
    std::string formatTimePretty(double time) const {
        static constexpr double thousand = 1000;
        static constexpr double million = thousand*thousand;
        static constexpr double billion = million*thousand;

        if(time < thousand) {
            return std::to_string(time) + "ns";
        }
        else if(time < million) {
            return std::to_string(time/thousand) + "us";
        }
        else if(time < billion) {
            return std::to_string(time/million) + "ms";
        }
        else {
            return std::to_string(time/billion) + "s";
        }
    }
    
    std::string getTimeString() const {
        return "Total Time: " + formatTimePretty(runTimesNanos) + "\tAverage Time: " + formatTimePretty(runTimesNanos / runCount);
    }
    
    friend std::ostream& operator<<(std::ostream& stream, const TimingDetails& timingDetails) {
        stream << timingDetails.name << ":\tRun Count: "<< timingDetails.runCount << "\t"<< timingDetails.getTimeString();
        return stream;
    }
};

class Timer
{
public:
    static Timer& getInstance()
    {
        static Timer    instance;
        return instance;
    }
    
    //Deleted copy functions
    Timer(Timer const&)           = delete;
    void operator=(Timer const&)  = delete;

private:
    Timer() {}
    
    void printTimingResults() {
#ifdef USE_TIMING
        std::cout<<"\nTIMING RESULTS\n";
        for(size_t i = 0; i < m_timingDetails.size(); i++) {
            std::cout<<m_timingDetails[i]<<"\n";
        }
        std::cout<<"\n";
#endif
    }
    
    void clearClaims() {
        for(size_t i = 0; i < m_timingDetails.size(); i++) {
            m_timingDetails[i].clearClaim();
        }
    }
    
    static size_t                   m_index;
    static bool                     m_enabled;
    std::vector<TimingDetails>      m_timingDetails;

public:
    inline void startTimer(size_t index) {
#ifdef USE_TIMING
        if(!Timer::isEnabled()) {
            return;
        }
        
        m_timingDetails[index].start = std::chrono::high_resolution_clock::now();
#endif
    }

    inline void stopTimer(size_t index) {
#ifdef USE_TIMING
        if(!Timer::isEnabled()) {
            return;
        }
        
        m_timingDetails[index].end = std::chrono::high_resolution_clock::now();
        m_timingDetails[index].runTimesNanos += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(m_timingDetails[index].end - m_timingDetails[index].start).count());
        m_timingDetails[index].runCount++;
#endif
    }
    
    inline void registerTimer(size_t& index, std::string name) {
#ifdef USE_TIMING
        if(!Timer::isEnabled()) {
            return;
        }
        
        if(index == -1) {
            index = m_index++;
            m_timingDetails.resize(m_index);
        }
                
        m_timingDetails[index].claim(name);
#endif
    }
    
    static void getTimingResultsAndClear() {
        getInstance().printTimingResults();
        getInstance().clearClaims();
    }
    
    static void enable() {
        m_enabled = true;
    }
    
    static void disable() {
        m_enabled = false;
    }
    
    static bool isEnabled() {
        return m_enabled;
    }
};
bool Timer::m_enabled = false;
size_t Timer::m_index = 0;




#endif /* utilities_h */
