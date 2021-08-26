//
//  Subdivision.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/27/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef Subdivision_Solver_h
#define Subdivision_Solver_h

#include "SolutionTracking/IntervalTracker.hpp"
#include "SolutionTracking/RootTracker.hpp"
#include "Utilities/utilities.hpp"
#include "Utilities/MultiPool.hpp"
#include "Utilities/ConcurrentStack.hpp"
#include "IntervalChecking/IntervalChecker.hpp"
#include "Approximation/ChebyshevApproximator.hpp"
#include "Solvers/LinearSolver.hpp"

template <int Rank>
class SubdivisionSolver
{
public:
    SubdivisionSolver(size_t _threadNum, const std::vector<Function::SharedFunctionPtr>& _functions, ConcurrentStack<SolveParameters>& _intervalsToRun, ObjectPool<SolveParameters>& _solveParametersPool, const SubdivisionParameters& _parameters, IntervalTracker& _intervalTracker, RootTracker& _rootTracker);
    ~SubdivisionSolver();

    void solve(SolveParameters* _parameters);
    
private:
    void subdivide(SolveParameters* _parameters, size_t _numGoodApproximations);
        
private:
    size_t                                                  m_threadNum;
    size_t                                                  m_rank;
    const std::vector<Function::SharedFunctionPtr>&         m_functions;
    ConcurrentStack<SolveParameters>&                       m_intervalsToRun;
    ObjectPool<SolveParameters>&                            m_solveParametersPool;
    SubdivisionParameters                                   m_subdivisionParameters;
    IntervalTracker&                                        m_intervalTracker;
    RootTracker&                                            m_rootTracker;
    std::vector<ChebyshevApproximation<Rank> >                 m_chebyshevApproximations;
    std::vector<std::unique_ptr<ChebyshevApproximator<Rank> > >m_chebyshevApproximators;
    IntervalChecker<Rank>                                      m_intervalChecker;
    LinearSolver<Rank>                                         m_linearSolver;
    std::vector<double>                                     m_minApproxTols;
    
    static size_t m_subdivisionSolverTimerIndex1;
    static size_t m_subdivisionSolverTimerIndex2;
    Timer& m_timer = Timer::getInstance();
};

template<int Rank>
size_t SubdivisionSolver<Rank>::m_subdivisionSolverTimerIndex1 = -1;

template<int Rank>
size_t SubdivisionSolver<Rank>::m_subdivisionSolverTimerIndex2 = -1;


#include "SubdivisionSolverND.ipp"

#endif /* Subdivision_Solver_h */
