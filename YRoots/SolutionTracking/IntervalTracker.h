//
//  IntervalData.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/6/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalData_h
#define IntervalData_h

#include "utilities.h"
#include <fstream>

struct IntervalResult {
    Interval    m_interval;
    SolveMethod m_howFound;
    
    IntervalResult() {}
    IntervalResult(Interval _interval, SolveMethod _howFound):
    m_interval(_interval),
    m_howFound(_howFound)
    {}
};

class IntervalTracker {
public:
    IntervalTracker(size_t _numThreads, bool _enabled): m_numThreads(_numThreads), m_enabled(_enabled)
    {
        m_intervalResults.resize(m_numThreads);
        m_outputFile = "intervals.txt";
    }
    
    void storeResult(size_t _threadNum, Interval& _interval, SolveMethod _howFound) {
        if (m_enabled) {
            m_intervalResults[_threadNum].emplace_back(_interval, _howFound);
        }
    }
    
    void printResults() {
        for(size_t i = 0; i < m_numThreads; i++) {
            std::cout<< "Thread " << i << " solved " << m_intervalResults[i].size() << " intervals\n";
        }
        return;
    }
    
    void logResults() {
        std::ofstream file;
        int precision = std::numeric_limits<double>::digits10 + 1;
        file.open (m_outputFile);
        for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
            file<<"Thread " << threadNum << " solved " << m_intervalResults[threadNum].size() << " intervals.\n";
            for(size_t intervalNum = 0; intervalNum < m_intervalResults[threadNum].size(); intervalNum++) {
                //Log the Lower Bound
                file<<"[";
                for(size_t i = 0; i < m_intervalResults[threadNum][intervalNum].m_interval.lowerBounds.size(); i++) {
                    file<<std::setprecision(precision)<<m_intervalResults[threadNum][intervalNum].m_interval.lowerBounds[i];
                    if(i + 1 < m_intervalResults[threadNum][intervalNum].m_interval.lowerBounds.size()){
                        file<<",";
                    }
                }
                file<<"],";
                //Log the Upper Bound
                file<<"[";
                for(size_t i = 0; i < m_intervalResults[threadNum][intervalNum].m_interval.upperBounds.size(); i++) {
                    file<<std::setprecision(precision)<<m_intervalResults[threadNum][intervalNum].m_interval.upperBounds[i];
                    if(i + 1 < m_intervalResults[threadNum][intervalNum].m_interval.upperBounds.size()){
                        file<<",";
                    }
                }
                file<<"]\t";
                //Log the reason
                file<<m_intervalResults[threadNum][intervalNum].m_howFound<<"\n";
            }
        }
        file.close();
    }

    
private:
    size_t                                      m_numThreads;
    bool                                        m_enabled;
    std::vector<std::vector<IntervalResult>>    m_intervalResults;
    std::string                                 m_outputFile;
};

#endif /* IntervalData_h */
