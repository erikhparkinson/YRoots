//
//  IntervalChecks.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/7/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalChecker_h
#define IntervalChecker_h

#include "../Approximation/ChebyshevApproximation.hpp"
#include "../SolutionTracking/IntervalTracker.h"
#include "../Utilities/MultiPool.h"
#include "../Utilities/ConcurrentStack.h"
#include "../Utilities/Timer.h"

enum EvalSign {
    //If we & to signs together as a bool, we get true if we can throw out an interval
    Zero = 0,
    Positive = 1,
    Negative = 2
};

template <Dimension D>
class IntervalChecker {
public:
    IntervalChecker(size_t _rank, IntervalTracker& _intervalTracker, size_t _threadNum, ConcurrentStack<SolveParameters>& _intervalsToRun, ObjectPool<SolveParameters>&        _solveParametersPool);
    
    bool runIntervalChecks(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval);
    void runSubintervalChecks(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations, SolveParameters* _currentParameters, size_t _numGoodApproximations);
    
protected:
    bool runConstantTermCheck(ChebyshevApproximation<D>& _approximation);
    void runQuadraticCheck(ChebyshevApproximation<D>& _approximation);
    double getBoundingInterval(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations);
    void pushIntervalToSolve(SolveParameters* _currentParameters, Interval& _newInterval);
    
    inline EvalSign getEvalSign(double _eval, double _error) {
        if(_eval > _error) {
            return EvalSign::Positive;
        }
        else if (_eval < -_error) {
            return EvalSign::Negative;
        }
        else {
            return EvalSign::Zero;
        }
    }

protected:
    size_t                  m_rank;
    IntervalTracker&        m_intervalTracker;
    static constexpr double m_randomIntervalDivider = 0.5139303900908738;
    std::vector<Interval>   m_scaledSubIntervals;
    std::vector<bool>       m_intervalMask;
    std::vector<bool>       m_throwOutMask;
    Interval                m_boundingInterval;
    Interval                m_tempInterval;
    
    //For BoundingIntervals
    std::vector<size_t>                             m_biggestArraySizeChecked; //Len Dim, one for each approximation
    std::vector<std::vector<std::vector<double>>>   m_reducedChebEvals; //Dim x (Dim-1) x length
                                                                        //For the x,y spot length is degree^y. Check degree against m_biggestArraySizeChecked each run.
    
    //Multithreading objects
    size_t                              m_threadNum;
    ConcurrentStack<SolveParameters>&   m_intervalsToRun;
    ObjectPool<SolveParameters>&        m_solveParametersPool;
    
    //For Timing
    static size_t           m_timerBoundingIntervalIndex;
    static size_t           m_timerQuadraticCheckIndex;
    Timer&                  m_timer = Timer::getInstance();
};

template<Dimension D>
size_t IntervalChecker<D>::m_timerBoundingIntervalIndex = -1;
template<Dimension D>
size_t IntervalChecker<D>::m_timerQuadraticCheckIndex = -1;

#include "BoundingIntervalUtilities.h"
#include "IntervalChecker1D.ipp"
#include "IntervalChecker2D.ipp"
#include "IntervalChecker3D.ipp"
#include "IntervalCheckerND.ipp"

#endif /* IntervalChecker_h */
