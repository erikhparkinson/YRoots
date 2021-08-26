//
//  ChebyshevApproximation1D.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/14/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ChebyshevApproximation1D_ipp
#define ChebyshevApproximation1D_ipp

template <>
void ChebyshevApproximation<1>::sumAbsValues() {
    if(!m_absValWasSummed){
        for(size_t i = 0; i < m_partialSideLength; i++) {
            m_sumAbsVal += std::abs(m_approximation[i]);
        }
        m_absValWasSummed = true;
    }
}

template <>
void ChebyshevApproximation<2>::sumAbsValues() {
    if(!m_absValWasSummed){
        for(size_t i = 0; i < m_partialSideLength; i++) {
            for(size_t j = 0; j < m_partialSideLength; j++) {
                m_sumAbsVal += std::abs(m_approximation[i + m_sideLength*j]);
            }
        }
        m_absValWasSummed = true;
    }
}

template <>
void ChebyshevApproximation<3>::sumAbsValues() {
    if(!m_absValWasSummed){
        for(size_t i = 0; i < m_partialSideLength; i++) {
            for(size_t j = 0; j < m_partialSideLength; j++) {
                for(size_t k = 0; k < m_partialSideLength; k++) {
                    m_sumAbsVal += std::abs(m_approximation[i + m_sideLength*(j + m_sideLength*k)]);
                }
            }
        }
        m_absValWasSummed = true;
    }
}

template <>
void ChebyshevApproximation<1>::setApproximation(size_t _degree, size_t _sideLength, double* _appoximation, double _infNorm, bool _signChange, double _approximationError)
{
    //Specialized to not set m_degreeSpots as they aren't used.
    m_partialSideLength = _degree+1;
    m_degree = _degree*m_rank;
    m_sideLength = _sideLength;
    m_approximation = _appoximation;
    m_infNorm = _infNorm;
    m_signChange = _signChange;
    m_approximationError = _approximationError;
    clear();
}

template <>
bool ChebyshevApproximation<1>::trimCoefficients(double _absApproxTol, double _relApproxTol, size_t _targetDegree) {
    //Continue while the approximation is good. This call updates goodDegree as well.
    while(isGoodApproximationSetDegree(_absApproxTol, _relApproxTol)) {
        if(m_degree <= _targetDegree) {
            return true;
        }
        if(m_absValWasSummed) {
            m_sumAbsVal -= std::abs(m_approximation[m_degree]);
        }
        m_approximationError += std::abs(m_approximation[m_degree]);
        m_approximation[m_degree] = 0;
        m_degree--;
    }
    return false;
}


#endif /* ChebyshevApproximation1D_ipp */
