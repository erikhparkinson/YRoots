//
//  LinearSolverND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef LinearSolverND_ipp
#define LinearSolverND_ipp

template <Dimension D>
LinearSolver<D>::LinearSolver(size_t _rank, IntervalTracker& _intervalTracker, RootTracker& _rootTracker) :
m_rank(_rank),
m_intervalTracker(_intervalTracker),
m_rootTracker(_rootTracker)
{
    m_linears = Eigen::MatrixXd::Zero(m_rank, m_rank);
    m_constants = Eigen::VectorXd::Zero(m_rank);
    m_result = Eigen::VectorXd::Zero(m_rank);
    m_root.resize(m_rank);
}

template <Dimension D>
void LinearSolver<D>::solve(std::vector<ChebyshevApproximation<D>>& _chebyshevApproximations, Interval& _interval, double _goodZerosTol) {
    //TODO: Implement this
}

#endif /* LinearSolverND_ipp */
