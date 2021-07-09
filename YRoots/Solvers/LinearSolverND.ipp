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
    //The i,j spot in the m_linears matrix is the jth linear term of the ith approximation
    //The ith spot in the m_constants vector is the 0th spot in the ith approximation.
    
    //We divide each equation by the inf norm of the function to make them similar sizes.
    //TODO: Is this the right way to handle this? Or should I just fix how colPivHouseholderQr handles it so it assumes it's full rank?
    
    for(size_t i = 0; i < m_rank; i++) { //Iterate through each approximation
        size_t spot = 1;
        size_t sideLength = _chebyshevApproximations[i].getSideLength();
        for(size_t j = 0; j < m_rank; j++) { //Iterate through the linear terms in the approximation
            m_linears(i,j) = _chebyshevApproximations[i].getArray()[spot] / _chebyshevApproximations[i].getInfNorm();
            spot *= sideLength;
        }
        //Set the constant term
        m_constants(i) = -_chebyshevApproximations[i].getArray()[0] / _chebyshevApproximations[i].getInfNorm();
    }

    m_result = m_linears.colPivHouseholderQr().solve(m_constants);
    for(size_t i = 0; i < m_rank; i++) {
        m_root[i] = m_result[i];
    }

    m_intervalTracker.storeResult(m_threadNum, _interval, SolveMethod::LinearSolve);
    m_rootTracker.storeRoot(m_threadNum, m_root, _interval, SolveMethod::LinearSolve, std::nan(""), _goodZerosTol);
}

#endif /* LinearSolverND_ipp */
