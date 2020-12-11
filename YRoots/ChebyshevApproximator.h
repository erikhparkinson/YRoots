//
//  ChebyshevApproximator.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/5/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ChebyshevApproximator_h
#define ChebyshevApproximator_h

#include "IntervalApproximator.h"
#include "ChebyshevApproximation.hpp"

template <Dimension D>
class ChebyshevApproximator
{
public:
    ChebyshevApproximator(size_t _rank, size_t _maxApproximationDegree, ChebyshevApproximation<D>& _approximation);
    ~ChebyshevApproximator();
    
    void approximate(const std::unique_ptr<Function>& _function, const Interval& _currentInterval, size_t _approximationDegree);
    
    bool hasSignChange() {
        return m_signChange;
    }
    
    ChebyshevApproximation<D>& getApproximation() {
        return m_approximation;
    }

private:
    void calculateApproximationError();
    
private:
    size_t                                  m_rank;
    size_t                                  m_maxApproximationDegree;
    
    double*                                 m_input;
    double*                                 m_output1;
    double*                                 m_output2;
    fftw_r2r_kind*                          m_kinds;
    double*                                 m_inputPartial;
    
    std::vector<std::unique_ptr<IntervalApproximator<D>>>    m_intervalApproximators;
    size_t                                  m_firstApproximator;
    size_t                                  m_secondApproximator;
    size_t                                  m_sideLength1;
    size_t                                  m_sideLength2;
    size_t                                  m_approxLength1;
    size_t                                  m_approxLength2;

    double                                  m_infNorm;
    bool                                    m_signChange;
    double                                  m_approximationError;
    
    ChebyshevApproximation<D>&              m_approximation;
};

    
#include "ChebyshevApproximatorND.ipp"

#endif /* ChebyshevApproximator_h */
