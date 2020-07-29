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
#include <tbb/concurrent_vector.h>

//THIS CLASS MUST BE THREAD SAFE

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
    IntervalTracker(bool _enabled): m_enabled(_enabled) {
        
    }
    
    void storeResult(Interval& _interval, SolveMethod _howFound) {
        if (m_enabled) {
            m_intervalResults.emplace_back(_interval, _howFound);
        }
    }
    
    void printResults() {
        std::cout<< m_intervalResults.size() << " Intervals Used:\n";
        return;
        for(size_t intervalNum = 0; intervalNum < m_intervalResults.size(); intervalNum++) {
            for(size_t i = 0; i < m_intervalResults[intervalNum].m_interval.lowerBounds.size(); i++) {
                std::cout<<m_intervalResults[intervalNum].m_interval.lowerBounds[i]<<",";
            }
            std::cout<<";\t";
            for(size_t i = 0; i < m_intervalResults[intervalNum].m_interval.upperBounds.size(); i++) {
                std::cout<<m_intervalResults[intervalNum].m_interval.upperBounds[i]<<",";
            }
            std::cout<<";\t";
            
            std::cout<<m_intervalResults[intervalNum].m_howFound<<"\n";
        }
    }
    
private:
    bool                                         m_enabled;
    tbb::concurrent_vector<IntervalResult>       m_intervalResults;
    
};

#endif /* IntervalData_h */
