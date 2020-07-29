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
void ChebyshevApproximation<Dimension::One>::sumAbsValues() {
    if(!m_absValWasSummed){
        for(size_t i = 0; i < m_degree+1; i++) {
            m_sumAbsVal += std::abs(m_approximation[i]);
        }
        m_absValWasSummed = true;
    }
}

template <>
void ChebyshevApproximation<Dimension::Two>::sumAbsValues() {
    if(!m_absValWasSummed){
        for(size_t i = 0; i < m_degree+1; i++) {
            for(size_t j = 0; j < m_degree+1; j++) {
                m_sumAbsVal += std::abs(m_approximation[i + m_sideLength*j]);
            }
        }
        m_absValWasSummed = true;
    }
}

template <>
void ChebyshevApproximation<Dimension::Three>::sumAbsValues() {
    if(!m_absValWasSummed){
        for(size_t i = 0; i < m_degree+1; i++) {
            for(size_t j = 0; j < m_degree+1; j++) {
                for(size_t k = 0; k < m_degree+1; k++) {
                    m_sumAbsVal += std::abs(m_approximation[i + m_sideLength*(j + m_sideLength*k)]);
                }
            }
        }
        m_absValWasSummed = true;
    }
}

template <>
bool ChebyshevApproximation<Dimension::One>::trimCoefficients(double _absApproxTol, double _relApproxTol, size_t _targetDegree) {
    //Continue while the approximation is good. This call updates goodDegree as well.
    while(isGoodApproximation(_absApproxTol, _relApproxTol)) {
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
