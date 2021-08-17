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
    double computeBoundingInterval(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations);

    const Interval& getBoundingInterval() {
        return m_boundingInterval;
    }
    
protected:
    double updateBoundingIntervalLinearErrorSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations);
    double updateBoundingIntervalLipshitzSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations);

    void preconditionPolynomials(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations);
    
    double computeLipshitzConstant(Eigen::VectorXd poly, size_t dim);
protected:
    size_t                  m_rank;
    Interval                m_boundingInterval;
            
    //For Bounding Intervals
    typename EigenTypes<Rank>::Matrix m_linearTerms;
    typename EigenTypes<Rank>::Vector m_constantTerms;
    typename EigenTypes<Rank>::Vector m_errorTerms;
    typename EigenTypes<Rank>::Vector m_errorTermsOfLinear;
    typename EigenTypes<Rank>::Matrix2Columns m_rightHandSideOfLinearSystemErrors;
    typename Eigen::ColPivHouseholderQR<typename EigenTypes<Rank>::Matrix> m_linearTermsQR;
    typename EigenTypes<Rank>::Matrix2Columns m_linearSystemWithErrorResult;
    
    //The Preconditioned Polynomials
    size_t m_preconditionPolysDegree;
    typename EigenTypes<Rank>::Vector m_preconditionedErrors;
    Eigen::MatrixXd     m_preconditionedPolys;
    
    //For Timing
    Timer&                  m_timer = Timer::getInstance();
};

#include "IntervalBounder1D.ipp"
#include "IntervalBounderND.ipp"

#endif /* IntervalBounder_h */
