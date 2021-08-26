//
//  LinearSolver.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef LinearSolver_h
#define LinearSolver_h

#include <complex>
#include "Eigen/Core"
#include "Eigen/Dense"

#include "Approximation/ChebyshevApproximation.hpp"
#include "Utilities/utilities.hpp"
#include "SolutionTracking/IntervalTracker.hpp"
#include "SolutionTracking/RootTracker.hpp"

template <int Rank>
class LinearSolver {
public:
    LinearSolver(size_t _threadNum, size_t _rank, IntervalTracker& _intervalTracker, RootTracker& _rootTracker);
    
    void solve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations, Interval& _interval, double _goodZerosTol);

private:
    size_t                  m_threadNum;
    size_t                  m_rank;
    IntervalTracker&        m_intervalTracker;
    RootTracker&            m_rootTracker;
    
    typename EigenTypes<Rank>::Matrix m_linearTerms;
    typename EigenTypes<Rank>::Vector m_constantTerms;
    typename EigenTypes<Rank>::Vector m_result;
    
    std::vector<std::complex<double> >     m_root;
};

#include "LinearSolver1D.ipp"
#include "LinearSolver2D.ipp"
#include "LinearSolverND.ipp"

#endif /* LinearSolver_h */
