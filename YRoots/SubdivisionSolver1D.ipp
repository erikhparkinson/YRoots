//
//  SubdivisionSolver1D.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/6/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef SubdivisionSolver1D_ipp
#define SubdivisionSolver1D_ipp

template <>
void SubdivisionSolver<Dimension::One>::solve(Interval _currentInterval, size_t currentLevel)
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


#endif /* SubdivisionSolver1D_ipp */
