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
SubdivisionSolver<D>::SubdivisionSolver(const std::vector<std::unique_ptr<FunctionInterface>>& _functions) :
m_functions(_functions),
m_rank(m_functions.size())
{
    for(size_t i = 0; i < m_rank; i++) {
        m_chebyshevApproximators.emplace_back(std::make_unique<ChebyshevApproximator<D>>(m_functions[i], m_subdivisionParameters.approximationDegree));
    }
}

template <Dimension D>
void SubdivisionSolver<D>::solve(Interval _currentInterval, size_t currentLevel)
{
    //If we are too deep, track that and return.
    if (currentLevel > m_subdivisionParameters.maxLevel)
    {
        m_intervalData.storeResult(_currentInterval, "Too Deep");
        return;
    }
    
    for(size_t funcNum = 0; funcNum < m_rank; funcNum++) {
        //Get a chebyshev approximation
        m_chebyshevApproximators[funcNum]->approximate(_currentInterval);
        if(!m_chebyshevApproximators[funcNum]->isGoodApproximation(m_subdivisionParameters)) {
            //Subdivide and continue
        }
        else if(!m_chebyshevApproximators[funcNum]->hasSignChange()) {
            //Run Checks
        }
    }
    
    //Trim the coeffs
        
}


#endif /* SubdivisionSolverND_ipp */
