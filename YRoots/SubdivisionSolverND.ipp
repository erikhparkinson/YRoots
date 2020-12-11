//
//  SubdivisionSolverND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/6/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef SubdivisionSolverND_ipp
#define SubdivisionSolverND_ipp

template <Dimension D>
SubdivisionSolver<D>::SubdivisionSolver(size_t _threadNum, const std::vector<std::unique_ptr<Function>>& _functions, ConcurrentStack<SolveParameters>& _intervalsToRun, ObjectPool<SolveParameters>& _solveParametersPool, SubdivisionParameters& _parameters, IntervalTracker& _intervalTracker, RootTracker& _rootTracker) :
m_threadNum(_threadNum),
m_rank(_functions.size()),
m_functions(_functions),
m_intervalsToRun(_intervalsToRun),
m_solveParametersPool(_solveParametersPool),
m_subdivisionParameters(_parameters),
m_intervalTracker(_intervalTracker),
m_rootTracker(_rootTracker),
m_intervalChecker(m_rank, m_intervalTracker, m_threadNum, m_intervalsToRun, m_solveParametersPool),
m_linearSolver(m_threadNum, m_rank, m_intervalTracker, m_rootTracker)
{
    for(size_t i = 0; i < m_functions.size(); i++) {
        //Create the Chebyshev Approximations
        m_chebyshevApproximations.emplace_back(m_rank);
    }

    for(size_t i = 0; i < m_functions.size(); i++) {
        //Create the Chebyshev Approximators
        m_chebyshevApproximators.emplace_back(std::make_unique<ChebyshevApproximator<D>>(m_rank, m_subdivisionParameters.approximationDegree, m_chebyshevApproximations[i]));
    }
}

template <Dimension D>
SubdivisionSolver<D>::~SubdivisionSolver(){

}

template <Dimension D>
void SubdivisionSolver<D>::subdivide(SolveParameters* _parameters, size_t _numGoodApproximations)
{
    m_intervalChecker.runSubintervalChecks(m_chebyshevApproximations, _parameters, _numGoodApproximations);
}

template <Dimension D>
void SubdivisionSolver<D>::solve(SolveParameters* _parameters)
{
    //Figure out how a speicific number is being solved
    /*if(_parameters->interval.lowerBounds[0] > 0.5176325837886175 || _parameters->interval.upperBounds[0] < 0.5176325837886175) {
        return;
    }
    if(_parameters->interval.lowerBounds[1] > 0.9659361387156291 || _parameters->interval.upperBounds[1] < 0.9659361387156291) {
        return;
    }*/

    //If we are too deep, track that and return.
    if (_parameters->currentLevel > m_subdivisionParameters.maxLevel)
    {
        m_intervalTracker.storeResult(m_threadNum, _parameters->interval, SolveMethod::TooDeep);
        return;
    }
    for(size_t funcNum = 0; funcNum < m_functions.size(); funcNum++) {
        //Get a chebyshev approximation
        m_chebyshevApproximators[funcNum]->approximate(m_functions[funcNum], _parameters->interval, _parameters->goodDegrees[funcNum]);
        if(!m_chebyshevApproximations[funcNum].isGoodApproximation(m_subdivisionParameters.absApproxTol, m_subdivisionParameters.relApproxTol)) {
            //Increase the goodDegree by 1 up to the max.
            _parameters->goodDegrees[funcNum] = std::min(_parameters->goodDegrees[funcNum]+1, m_subdivisionParameters.approximationDegree);
            //Subdivide
            return subdivide(_parameters, funcNum);
        }
        else if(!m_chebyshevApproximations[funcNum].hasSignChange()) {
            //Run Checks
            if(!m_intervalChecker.runIntervalChecks(m_chebyshevApproximations[funcNum], _parameters->interval)){
                return;
            }
        }
        //Update the good degree
        _parameters->goodDegrees[funcNum] = std::min(_parameters->goodDegrees[funcNum], m_chebyshevApproximations[funcNum].getGoodDegree());
    }
    
    //Trim the coeffs
    bool goodApproximations = true;
    for(size_t funcNum = 0; funcNum < m_functions.size(); funcNum++) {
        //Trim
        goodApproximations &= m_chebyshevApproximations[funcNum].trimCoefficients(m_subdivisionParameters.absApproxTol, m_subdivisionParameters.relApproxTol, m_subdivisionParameters.targetDegree);
        //Get the good degree
        _parameters->goodDegrees[funcNum] = std::min(_parameters->goodDegrees[funcNum], m_chebyshevApproximations[funcNum].getGoodDegree());
    }
    if(!goodApproximations) {
        return subdivide(_parameters, m_functions.size());
    }

    //Check if everything is linear and get goodZerosTol
    bool isLinear = true;
    double goodZerosTol = 0;
    for(size_t funcNum = 0; funcNum < m_functions.size(); funcNum++) {
        if(!m_chebyshevApproximations[funcNum].isLinear()) {
            isLinear = false;
        }
        goodZerosTol += m_chebyshevApproximations[funcNum].getApproximationError();
    }
    goodZerosTol = std::max(m_subdivisionParameters.minGoodZerosTol, goodZerosTol*m_subdivisionParameters.goodZerosFactor);
    //Solve
    if(isLinear) {
        return m_linearSolver.solve(m_chebyshevApproximations, _parameters->interval, goodZerosTol);
    }
    else {
        //TODO: Sove using spectral methods
        return subdivide(_parameters, m_functions.size());
    }
}


#endif /* SubdivisionSolverND_ipp */
