//
//  LinearSolver1D.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef LinearSolver1D_ipp
#define LinearSolver1D_ipp

template <>
void LinearSolver<1>::solve(std::vector<ChebyshevApproximation<1> >& _chebyshevApproximations, Interval& _interval, double _goodZerosTol) {
    m_root[0] = -_chebyshevApproximations[0].getArray()[0]/_chebyshevApproximations[0].getArray()[1];
    const bool hasRoot = m_rootTracker.storeRoot(m_threadNum, m_root, _interval, SolveMethod::LinearSolve, std::nan(""), _goodZerosTol);
    m_intervalTracker.storeResult(m_threadNum, _interval, SolveMethod::LinearSolve, hasRoot);
}

#endif /* LinearSolver1D_ipp */
