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

template <Dimension D>
class ChebyshevApproximator
{
public:
    ChebyshevApproximator(const std::unique_ptr<FunctionInterface>& _function, size_t _approximationDegree);
    ~ChebyshevApproximator();
    
    void approximate(const Interval& _currentInterval);
    bool isGoodApproximation(SubdivisionParameters& _subdivisionParameters);
    
    bool hasSignChange() {
        return m_signChange;
    }
    
    double* getApproximation() {
        return m_intervalApproximator1.getOutput();
    }

private:
    void calculateApproximationError();
    
private:
    IntervalApproximator<D> m_intervalApproximator1;
    IntervalApproximator<D> m_intervalApproximator2;
    size_t                  m_rank;
    size_t                  m_approximationDegree;
    size_t                  m_sideLength1;
    size_t                  m_sideLength2;
    double                  m_infNorm;
    bool                    m_signChange;
    double                  m_approximationError;
};

    
#include "ChebyshevApproximator1D.ipp"
#include "ChebyshevApproximatorND.ipp"

#endif /* ChebyshevApproximator_h */
