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
#include "MultiPool.h"
#include "ConcurrentStack.h"

template <Dimension D>
class IntervalChecker {
public:
    IntervalChecker(size_t _rank, IntervalTracker& _intervalTracker, size_t _threadNum, ConcurrentStack<SolveParameters>& _intervalsToRun, ObjectPool<SolveParameters>&        _solveParametersPool);
    
    bool runIntervalChecks(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval);
    void runSubintervalChecks(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations, SolveParameters* _currentParameters, size_t _numGoodApproximations);
    
private:
    bool runConstantTermCheck(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval);
    bool runQuadraticCheck(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval);

private:
    size_t                  m_rank;
    IntervalTracker&        m_intervalTracker;
    double                  m_randomIntervalDivider;
    std::vector<Interval>   m_scaledSubIntervals;
    
    //Multithreading objects
    size_t                              m_threadNum;
    ConcurrentStack<SolveParameters>&   m_intervalsToRun;
    ObjectPool<SolveParameters>&        m_solveParametersPool;

    
};

#include "IntervalChecker1D.ipp"
#include "IntervalChecker2D.ipp"
#include "IntervalChecker3D.ipp"
#include "IntervalCheckerND.ipp"

#endif /* IntervalChecker_h */
