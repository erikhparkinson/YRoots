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
IntervalChecker<D>::IntervalChecker(IntervalTracker& _intervalTracker, size_t _rank, tbb::strict_ppl::concurrent_queue<SolveParameters>& _intervalsToRun):
m_rank(_rank),
m_intervalTracker(_intervalTracker),
m_randomIntervalDivider(0.5139303900908738),
m_intervalsToRun(_intervalsToRun)
{
    //Create the scaled subintervals
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
    }
}

template <Dimension D>
bool IntervalChecker<D>::runIntervalChecks(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval)
{
    //Returns true if we need to keep the interval
    if(!runConstantTermCheck(_approximation, _currentInterval)) {
        m_intervalTracker.storeResult(_currentInterval, SolveMethod::ConstantTermCheck);
        return false;
    }
    else {
        return true;
    }
}

template <Dimension D>
void IntervalChecker<D>::runSubintervalChecks(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations, SolveParameters& _currentParameters, size_t _numGoodApproximations)
{
    for(size_t intervalNum = 0; intervalNum < m_scaledSubIntervals.size(); intervalNum++) {
        SolveParameters _nextParameters = _currentParameters;
        _nextParameters.currentLevel++;
        
        //Transform this interval onto m_scaledSubIntervals[intervalNum]
        for(size_t i = 0; i < m_rank; i++) {
            
            //TODO: Make some function for projections. Make this more efficient.
            double temp1 = _currentParameters.interval.upperBounds[i] - _currentParameters.interval.lowerBounds[i];
            double temp2 = _currentParameters.interval.upperBounds[i] + _currentParameters.interval.lowerBounds[i];

            _nextParameters.interval.lowerBounds[i] = ((temp1 * m_scaledSubIntervals[intervalNum].lowerBounds[i]) + temp2) /2.0;
            _nextParameters.interval.upperBounds[i] = ((temp1 * m_scaledSubIntervals[intervalNum].upperBounds[i]) + temp2) /2.0;
        }
        
        //TODO: Run the checks here
        //It's probably the most efficient to run the checks on all the intervals at the same time to avoid
        //repeating calculations. Either that, or find a way to store calculations inside of the ChebyshevApproximation.
        
        m_intervalsToRun.push(_nextParameters);
    }
    
    
    //TODO: Write this function!
    //TODO: Have this class hold a reference to interval data and store the result
    /*size_t currSpot = 0;
    while (currSpot < _subIntervals.size()) {
        if(!runQuadraticCheck(_approximation, _subIntervals[currSpot])) {
            m_intervalData.storeResult(_subIntervals[currSpot], "QuadraticCheck");
            _subIntervals.erase(_subIntervals.begin() + currSpot);
        }
        else {
            currSpot++;
        }
    }*/
    
    
    
    //for(size_t i = 0; i < m_scaledSubIntervals.size(); i++) {
        //Construct the new interval
        //Interval newInterval;
        /*for(size_t j = 0; j < m_rank; j++) {
            newInterval.lowerBounds.push_back(_currentInterval.lowerBounds[j]*m_scaledSubIntervals[i].lowerBounds[j]);
            newInterval.upperBounds.push_back(_currentInterval.upperBounds[j]*m_scaledSubIntervals[i].upperBounds[j]);
        }*/
        //Run the subinterval tests
        //if(!m_intervalChecker.runSubintervalChecks(m_chebyshevApproximators[funcNum]->getApproximation(), _currentInterval)) {
        //    solve(_currentInterval, _currentLevel+1);
        //}
    //}
    
    //Run subinterval the tests on it.
    //Subdivide if needed.

}

template <Dimension D>
bool IntervalChecker<D>::runConstantTermCheck(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval)
{
    //Returns true if we need to keep the interval
    _approximation.sumAbsValues();
    return _approximation.getSumAbsVal() + _approximation.getApproximationError() > std::abs(2*_approximation.getArray()[0]);
}

template <Dimension D>
bool IntervalChecker<D>::runQuadraticCheck(ChebyshevApproximation<D>& _approximation, Interval& _currentInterval)
{
    //Returns true if we need to keep the interval
    //TODO: Write this function!
    return true;
}



#endif /* IntervalCheckerND_h */
