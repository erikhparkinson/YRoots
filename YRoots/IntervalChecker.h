//
//  IntervalChecks.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/7/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalChecker_h
#define IntervalChecker_h

#include "ChebyshevApproximation.hpp"
#include "IntervalTracker.h"

template <Dimension D>
class IntervalChecker {
public:
    IntervalChecker(IntervalTracker& _intervalTracker, size_t _rank, tbb::strict_ppl::concurrent_queue<SolveParameters>& _intervalsToRun);
    
    bool runIntervalChecks(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval);
    void runSubintervalChecks(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations, SolveParameters& _currentParameters, size_t _numGoodApproximations);
    
private:
    bool runConstantTermCheck(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval);
    bool runQuadraticCheck(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval);

private:
    size_t                  m_rank;
    IntervalTracker&        m_intervalTracker;
    double                  m_randomIntervalDivider;
    std::vector<Interval>   m_scaledSubIntervals;
    
    tbb::strict_ppl::concurrent_queue<SolveParameters>&     m_intervalsToRun;
};

#include "IntervalChecker1D.ipp"
#include "IntervalChecker2D.ipp"
#include "IntervalChecker3D.ipp"
#include "IntervalCheckerND.ipp"

#endif /* IntervalChecker_h */
