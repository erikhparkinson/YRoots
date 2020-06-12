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
    ChebyshevApproximator(const std::unique_ptr<FunctionInterface>& _function, size_t approximationDegree);
    ~ChebyshevApproximator();
    
    void approximate(const Interval& _currentInterval);
    
private:
    IntervalApproximator<D> m_intervalApproximator1;
    IntervalApproximator<D> m_intervalApproximator2;
};

    
#include "ChebyshevApproximator1D.ipp"
#include "ChebyshevApproximatorND.ipp"

#endif /* ChebyshevApproximator_h */
