//
//  IntervalCheckerND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/11/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalCheckerND_h
#define IntervalCheckerND_h

template <Dimension D>
IntervalChecker<D>::IntervalChecker(size_t _rank, IntervalTracker& _intervalTracker, size_t _threadNum, ConcurrentStack<SolveParameters>& _intervalsToRun, ObjectPool<SolveParameters>& _solveParametersPool):
m_rank(_rank),
m_intervalTracker(_intervalTracker),
m_threadNum(_threadNum),
m_intervalsToRun(_intervalsToRun),
m_solveParametersPool(_solveParametersPool)
{
    //Create the scaled subintervals
    //[-,-,...,-,-]
    //[-,-,...,-,+]
    //[-,-,...,+,-]
    //[-,-,...,+,+]
    //............
    //[+,+,...,+,+]
    double midPoint = 2*m_randomIntervalDivider - 1;
    size_t numIntervals = (1 << m_rank);
    for(size_t i = 0 ; i < numIntervals; i++) {
        Interval tempInterval;
        for(size_t j = m_rank-1; j < m_rank; j--) {
            if((i>>j)%2) {
                tempInterval.lowerBounds.push_back(midPoint);
                tempInterval.upperBounds.push_back(1.0);
            }
            else{
                tempInterval.lowerBounds.push_back(-1.0);
                tempInterval.upperBounds.push_back(midPoint);
            }
        }
        m_scaledSubIntervals.push_back(tempInterval);
        m_intervalMask.push_back(false);
        m_throwOutMask.push_back(false);
    }
    //Initialize the boudning interval
    for(size_t i = 0; i < m_rank; i++) {
        m_boundingInterval.lowerBounds.push_back(0.0);
        m_boundingInterval.upperBounds.push_back(0.0);
        m_tempInterval.lowerBounds.push_back(0.0);
        m_tempInterval.upperBounds.push_back(0.0);
    }
    
    m_timer.registerTimer(m_timerBoundingIntervalIndex, "Bounding Interval");
    m_timer.registerTimer(m_timerQuadraticCheckIndex, "Quadratic Check");
}

template <Dimension D>
bool IntervalChecker<D>::runIntervalChecks(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval)
{
    //Returns true if we need to keep the interval
    if(!runConstantTermCheck(_approximation)) {
        m_intervalTracker.storeResult(m_threadNum, _currentInterval, SolveMethod::ConstantTermCheck);
        return false;
    }
    else {
        return true;
    }
}

template <Dimension D>
void IntervalChecker<D>::runSubintervalChecks(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations, SolveParameters* _currentParameters, size_t _numGoodApproximations)
{
    //TODO: Do I need this or should I just use the size in the Interval???
    m_timer.startTimer(m_timerBoundingIntervalIndex);
    double boundingIntervalSize = getBoundingInterval(_chebyshevApproximations);
    m_timer.stopTimer(m_timerBoundingIntervalIndex);
    
    if(boundingIntervalSize == 0.0) {
        //Track the interval and return
        m_intervalTracker.storeResult(m_threadNum, _currentParameters->interval, SolveMethod::BoundingInterval);
        return;
    }
    else if(boundingIntervalSize < 1.0) {
        //Track it and push it
        m_intervalTracker.storeResult(m_threadNum, _currentParameters->interval, SolveMethod::BoundingInterval, boundingIntervalSize);
        pushIntervalToSolve(_currentParameters, m_boundingInterval);
        return;
    }
    
    //Reset the mask
    for(size_t i = 0; i < m_scaledSubIntervals.size(); i++) {
        m_intervalMask[i] = false;
    }
    
    //Run the quad checks
    for(size_t i = 0; i < _numGoodApproximations; i++) {
        m_timer.startTimer(m_timerQuadraticCheckIndex);
        runQuadraticCheck(_chebyshevApproximations[i]);
        m_timer.stopTimer(m_timerQuadraticCheckIndex);
    }
    
    //Check how many intervals we are running
    size_t intervalsToKeep = 0;
    for(size_t i = 0; i < m_intervalMask.size(); i++) {
        intervalsToKeep += m_intervalMask[i];
    }
    
    //Use the bounding interval if it's better
    if(boundingIntervalSize < intervalsToKeep && boundingIntervalSize < 1.5) {
        //Track it and push it
        m_intervalTracker.storeResult(m_threadNum, _currentParameters->interval, SolveMethod::BoundingInterval, boundingIntervalSize);
        pushIntervalToSolve(_currentParameters, m_boundingInterval);
        return;
    }

    //Subdivide the intervals left over from the checks
    for(size_t intervalNum = 0; intervalNum < m_scaledSubIntervals.size(); intervalNum++) {
        if(m_intervalMask[intervalNum]) {
            //Get the projected interval and store it
            projectInterval(m_tempInterval, _currentParameters->interval, m_scaledSubIntervals[intervalNum]);
            m_intervalTracker.storeResult(m_threadNum, m_tempInterval, SolveMethod::QuadrticCheck);
        }
        else {
            //Push it to subdivision
            pushIntervalToSolve(_currentParameters, m_scaledSubIntervals[intervalNum]);
        }
    }
}

template <Dimension D>
bool IntervalChecker<D>::runConstantTermCheck(ChebyshevApproximation<D>& _approximation)
{
    //Returns true if we need to keep the interval
    _approximation.sumAbsValues();
    return _approximation.getSumAbsVal() + _approximation.getApproximationError() > std::abs(2*_approximation.getArray()[0]);
}

template <Dimension D>
void IntervalChecker<D>::runQuadraticCheck(ChebyshevApproximation<D>& _approximation)
{
    //TODO: Updates m_intervalMask to false for each interval we should remove
}

template <Dimension D>
double IntervalChecker<D>::getBoundingInterval(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations) {
    //TODO: Updates m_boundingInterval. Returns the size of the bounding interval
    return std::numeric_limits<double>::max();
}

template <Dimension D>
void IntervalChecker<D>::pushIntervalToSolve(SolveParameters* _currentParameters, Interval& _newInterval) {
    SolveParameters* _nextParameters = m_solveParametersPool.pop();
    //Get the projected interval
    projectInterval(_nextParameters->interval, _currentParameters->interval, _newInterval);

    //Populate the other details
    _nextParameters->currentLevel = _currentParameters->currentLevel+1;
    for(size_t i = 0; i < _nextParameters->goodDegrees.size(); i++) {
        _nextParameters->goodDegrees[i] = _currentParameters->goodDegrees[i];
    }
    
    //Push the new interval
    m_intervalsToRun.push(m_threadNum, _nextParameters);
}


#endif /* IntervalCheckerND_h */
