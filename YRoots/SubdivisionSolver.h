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
#include "IntervalTracker.h"
#include "RootTracker.h"
#include "ChebyshevApproximator.h"
#include "LinearSolver.h"

template <Dimension D>
class SubdivisionSolver
{
public:
    SubdivisionSolver(const std::vector<std::unique_ptr<FunctionInterface>>& _functions, tbb::strict_ppl::concurrent_queue<SolveParameters>& _intervalsToRun, SubdivisionParameters _parameters, IntervalTracker& _intervalTracker, RootTracker& _rootTracker);
    ~SubdivisionSolver();

    void solve(SolveParameters& _parameters);
    
private:
    void subdivide(SolveParameters& _parameters, size_t _numGoodApproximations);
        
private:
    size_t                                                  m_rank;
    const std::vector<std::unique_ptr<FunctionInterface>>&  m_functions;
    tbb::strict_ppl::concurrent_queue<SolveParameters>&     m_intervalsToRun;
    SubdivisionParameters                                   m_subdivisionParameters;
    IntervalTracker&                                        m_intervalTracker;
    RootTracker&                                            m_rootTracker;
    std::vector<ChebyshevApproximation<D>>                  m_chebyshevApproximations;
    std::vector<std::unique_ptr<ChebyshevApproximator<D>>>  m_chebyshevApproximators;
    IntervalChecker<D>                                      m_intervalChecker;
    LinearSolver<D>                                         m_linearSolver;
};

#include "SubdivisionSolverND.ipp"

#endif /* Subdivision_Solver_h */
