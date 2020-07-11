//
//  Subdivision.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/27/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef Subdivision_Solver_h
#define Subdivision_Solver_h

#include "utilities.h"
#include "IntervalChecker.h"
#include "IntervalData.h"
#include "ChebyshevApproximator.h"

template <Dimension D>
class SubdivisionSolver
{
public:
    SubdivisionSolver(const std::vector<std::unique_ptr<FunctionInterface>>& _functions);
    
    void solve(Interval _currentInterval, size_t currentLevel);
    
private:
    const std::vector<std::unique_ptr<FunctionInterface>>&  m_functions;
    SubdivisionParameters                                   m_subdivisionParameters;
    IntervalData                                            m_intervalData;
    std::vector<std::unique_ptr<ChebyshevApproximator<D>>>  m_chebyshevApproximators;
    size_t                                                  m_rank;
    
    
    //TODO: Create Interval Data Class
    //TODO: Create objects to handle root tracking
    //TODO: Pass in the tolerances
    //TODO: Make a class for FFTs
    //
    
    
};

#include "SubdivisionSolver1D.ipp"
#include "SubdivisionSolverND.ipp"

#endif /* Subdivision_Solver_h */
