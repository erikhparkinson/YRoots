//
//  ChebyshevApproximation.hpp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/9/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ChebyshevApproximation_hpp
#define ChebyshevApproximation_hpp

#include <stdio.h>

class ChebyshevApproximation
{
public:
    ChebyshevApproximation() :
    m_rank(0),
    m_degree(0),
    m_sideLength(0),
    m_approxLength(0),
    m_appoximation(nullptr),
    m_infNorm(0),
    m_signChange(false),
    m_approximationError(0)
    {
        
    }
    
    void setApproximation(size_t _rank, size_t _degree, size_t _sideLength, size_t _approxLength, double* _appoximation, double _infNorm, bool _signChange, double _approximationError)
    {
        m_rank = _rank;
        m_degree = _degree;
        m_sideLength = _sideLength;
        m_approxLength = _approxLength;
        m_appoximation = _appoximation;
        m_infNorm = _infNorm;
        m_signChange = _signChange;
        m_approximationError = _approximationError;
    }
    
    double* getArray() {
        return m_appoximation;
    }
    
private:
    size_t      m_rank;
    size_t      m_degree;
    size_t      m_sideLength;
    size_t      m_approxLength;
    double*     m_appoximation;
    
    double m_infNorm;
    bool m_signChange;
    double m_approximationError;
};

#endif /* ChebyshevApproximation_hpp */
