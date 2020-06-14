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
m_functions(_functions)
{
    for(size_t i = 0; i<m_functions.size(); i++) {
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
    
    //Get a chebyshev approximation
    std::cout<<m_functions[0]->to_string()<<"\n";
    m_chebyshevApproximators[0]->approximate(_currentInterval);
}


#endif /* SubdivisionSolverND_ipp */
