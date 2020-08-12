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
    
private:
    size_t                                      m_numThreads;
    bool                                        m_enabled;
    std::vector<std::vector<IntervalResult>>    m_intervalResults;
    
};

#endif /* IntervalData_h */
