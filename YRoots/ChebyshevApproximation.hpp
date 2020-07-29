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

template <Dimension D>
class ChebyshevApproximation
{
public:
    ChebyshevApproximation();
    
    void setApproximation(size_t _rank, size_t _degree, size_t _sideLength, double* _appoximation, double _infNorm, bool _signChange, double _approximationError);
    
    void clear();
    
    inline bool isGoodApproximation(double absApproxTol, double relApproxTol);
        
    bool trimCoefficients(double _absApproxTol, double _relApproxTol, size_t _targetDegree);
    
    void sumAbsValues();
    
    double* getArray();
    
    bool isLinear() {
        return m_degree == 1;
    }

    double getSumAbsVal() {
        return m_sumAbsVal;
    }
    
    bool hasSignChange() {
        return m_signChange;
    }
    
    double getApproximationError() {
        return m_approximationError;
    }
    
    size_t getGoodDegree() {
        return m_goodDegree;
    }
    
private:
    size_t      m_rank;
    size_t      m_degree;
    size_t      m_sideLength;
    double*     m_approximation;
    
    double      m_infNorm;
    bool        m_signChange;
    double      m_approximationError;
    
    bool        m_absValWasSummed;
    double      m_sumAbsVal;
    size_t      m_goodDegree;
};

#include "ChebyshevApproximation1D.ipp"
#include "ChebyshevApproximationND.ipp"

#endif /* ChebyshevApproximation_hpp */
