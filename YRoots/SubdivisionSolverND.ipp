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
SubdivisionSolver<D>::SubdivisionSolver(const std::vector<std::unique_ptr<FunctionInterface>>& _functions, tbb::strict_ppl::concurrent_queue<SolveParameters>& _intervalsToRun) :
m_rank(_functions.size()),
m_functions(_functions),
m_intervalsToRun(_intervalsToRun),
m_intervalChecker(m_intervalData, m_rank, m_intervalsToRun)
{
    for(size_t i = 0; i < m_rank; i++) {
        //Create the Chebyshev Approximations
        m_chebyshevApproximations.emplace_back();

        //Create the Chebyshev Approximators
        m_chebyshevApproximators.emplace_back(std::make_unique<ChebyshevApproximator<D>>(m_rank, m_subdivisionParameters.approximationDegree, m_chebyshevApproximations[i]));
    }
}

template <Dimension D>
SubdivisionSolver<D>::~SubdivisionSolver(){

}

template <Dimension D>
void SubdivisionSolver<D>::subdivide(SolveParameters& _parameters, size_t _numGoodApproximations)
{
    m_intervalChecker.runSubintervalChecks(m_chebyshevApproximations, _parameters, _numGoodApproximations);
}

template <Dimension D>
void SubdivisionSolver<D>::solve(SolveParameters& _parameters)
{
    //If we are too deep, track that and return.
    if (_parameters.currentLevel > m_subdivisionParameters.maxLevel)
    {
        m_intervalData.storeResult(_parameters.interval, "Too Deep");
        return;
    }
    
    for(size_t funcNum = 0; funcNum < m_rank; funcNum++) {
        //Get a chebyshev approximation
        m_chebyshevApproximators[funcNum]->approximate(m_functions[funcNum], _parameters.interval, m_subdivisionParameters.approximationDegree);
        if(!m_chebyshevApproximators[funcNum]->isGoodApproximation(m_subdivisionParameters)) {
            //Subdivide and continue
            subdivide(_parameters, funcNum);
        }
        else if(!m_chebyshevApproximators[funcNum]->hasSignChange()) {
            //Run Checks
            if(m_intervalChecker.runIntervalChecks(m_chebyshevApproximators[funcNum]->getApproximation(), _parameters.interval)){
                return;
            }
        }
    }
    
    //Trim the coeffs
    
    
    //If trim coeffs introduces too much error than subdivide
    
    
    //If any degree is too big then subdivide
    
    
    //Check if everything is linear
    
    
    //Sove using spectral methods
}


#endif /* SubdivisionSolverND_ipp */
