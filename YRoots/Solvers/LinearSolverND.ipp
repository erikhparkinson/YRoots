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
LinearSolver<D>::LinearSolver(size_t _threadNum, size_t _rank, IntervalTracker& _intervalTracker, RootTracker& _rootTracker) :
m_threadNum(_threadNum),
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
    //The j,i spot in the m_linears matrix is the ith linear term of the jth approximation
    //The ith spot in the m_constants vector is the 0th spot in the ith approximation.
    
    size_t spot = 1;
    size_t sideLength = _chebyshevApproximations[0].getSideLength();
    for(size_t i = 0; i < m_rank; i++) {
        //Set all the ith terms for each approximaitons.
        for(size_t j = 0; j < m_rank; j++) {
            m_linears(j,i) = _chebyshevApproximations[j].getArray()[spot];
        }
        spot *= sideLength;
        //Set the constant terms
        m_constants(i) = _chebyshevApproximations[i].getArray()[0];
    }
    
    m_result = m_linears.colPivHouseholderQr().solve(m_constants);
    for(size_t i = 0; i < m_rank; i++) {
        m_root[i] = m_result[i];
    }
    
    m_intervalTracker.storeResult(m_threadNum, _interval, SolveMethod::LinearSolve);
    m_rootTracker.storeRoot(m_threadNum, m_root, _interval, SolveMethod::LinearSolve, std::nan(""), _goodZerosTol);
}

#endif /* LinearSolverND_ipp */
