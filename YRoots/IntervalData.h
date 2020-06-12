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
    std::string m_howFound;
    
    IntervalResult() {}
    IntervalResult(Interval _interval, std::string _howFound):
    m_interval(_interval),
    m_howFound(_howFound)
    {}
};

class IntervalData {
public:
    IntervalData() {
        
    }
    
    void storeResult(Interval _interval, std::string _howFound) {
        m_intervalResults.emplace_back(_interval, _howFound);
    }
    
private:
    std::vector<IntervalResult>       m_intervalResults;
    
};

#endif /* IntervalData_h */
