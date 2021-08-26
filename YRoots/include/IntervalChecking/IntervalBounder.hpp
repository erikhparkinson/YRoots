//
//  IntervalBounder.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 8/14/21.
//  Copyright © 2021 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalBounder_h
#define IntervalBounder_h

#include "Approximation/ChebyshevApproximation.hpp"
#include "Utilities/Timer.hpp"
#include "Utilities/utilities.hpp"
#include <Eigen/Dense>

#include "IntervalChecking/BoundingIntervalUtilities.hpp"

template <int Rank>
class IntervalBounder {
public:
    IntervalBounder(size_t _rank);
    double computeBoundingInterval(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations, const std::vector<bool>& _allowedToReduceDim);

    const Interval& getBoundingInterval() {
        return m_boundingInterval;
    }
    
protected:
    double updateBoundingIntervalLinearErrorSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations, const std::vector<bool>& _allowedToReduceDim);
    double updateBoundingIntervalLipshitzSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations, const std::vector<bool>& _allowedToReduceDim);

    void preconditionPolynomials(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations);
    
    double computeLipshitzConstant(const Eigen::VectorXd& poly, size_t dim);
    void chebValReduce(const Eigen::VectorXd& poly, Eigen::VectorXd& result, size_t dim, double value);
    double getLiphsitzBoundIncreaseND(const Eigen::VectorXd& poly, double error, double lipshitzConstant, const Interval& boundingInterval, size_t dim);
    double getExtremeAbsVal(const Eigen::VectorXd& poly, const Interval& boundingInterval, size_t dim);
protected:
    size_t                  m_rank;
    Interval                m_boundingInterval;
            
    //For Bounding Intervals
    typename EigenTypes<Rank>::Matrix m_linearTerms;
    typename EigenTypes<Rank>::Vector m_constantTerms;
    typename EigenTypes<Rank>::Vector m_errorTerms;
    typename EigenTypes<Rank>::Vector m_errorTermsOfLinear;
    typename EigenTypes<Rank>::MatrixPowerColumns m_rightHandSideOfLinearSystemErrors;
    typename Eigen::ColPivHouseholderQR<typename EigenTypes<Rank>::Matrix> m_linearTermsQR;
    typename EigenTypes<Rank>::MatrixPowerColumns m_linearSystemWithErrorResult;
    typename EigenTypes<Rank>::Matrix m_linearTermsInverse;

    //The Preconditioned Polynomials
    size_t m_preconditionPolysDegree;
    typename EigenTypes<Rank>::Vector m_preconditionedErrors;
    Eigen::MatrixXd     m_preconditionedPolys;
    
    //For Timing
    static size_t           m_timerLinearErrorSolveIndex;
    static size_t           m_timerPreconditionPolysIndex;
    static size_t           m_timerLipshitzSolveIndex;
    static size_t           m_timerChebValReduceIndex;
    Timer&                  m_timer = Timer::getInstance();
};

template<int Rank>
size_t IntervalBounder<Rank>::m_timerLinearErrorSolveIndex = -1;
template<int Rank>
size_t IntervalBounder<Rank>::m_timerPreconditionPolysIndex = -1;
template<int Rank>
size_t IntervalBounder<Rank>::m_timerLipshitzSolveIndex = -1;
template<int Rank>
size_t IntervalBounder<Rank>::m_timerChebValReduceIndex = -1;

#include "IntervalBounder1D.ipp"
#include "IntervalBounderND.ipp"

#endif /* IntervalBounder_h */
