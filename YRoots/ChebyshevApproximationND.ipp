//
//  ChebyshevApproximationND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/14/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ChebyshevApproximationND_ipp
#define ChebyshevApproximationND_ipp

template <Dimension D>
ChebyshevApproximation<D>::ChebyshevApproximation():
m_rank(0),
m_degree(0),
m_sideLength(0),
m_approximation(nullptr),
m_infNorm(0),
m_signChange(false),
m_approximationError(0),
m_absValWasSummed(false),
m_sumAbsVal(0)
{
    
}

template <Dimension D>
void ChebyshevApproximation<D>::setApproximation(size_t _rank, size_t _degree, size_t _sideLength, double* _appoximation, double _infNorm, bool _signChange, double _approximationError)
{
    m_rank = _rank;
    m_degree = _degree;
    m_sideLength = _sideLength;
    m_approximation = _appoximation;
    m_infNorm = _infNorm;
    m_signChange = _signChange;
    m_approximationError = _approximationError;
}

template <Dimension D>
double* ChebyshevApproximation<D>::getArray() {
    return m_approximation;
}

template <Dimension D>
void ChebyshevApproximation<D>::sumAbsValues() {
    if(!m_absValWasSummed){
        //Set up the needed variables
        std::vector<size_t> inputSpot(m_rank,0);
        std::vector<size_t> multipliers(m_rank, 1);
        for(size_t i = 1; i < m_rank; i++) {
            multipliers[i] = m_sideLength*multipliers[i-1];
        }
        
        //Iterate through all the combinations
        size_t spotToInc = 0;
        m_sumAbsVal += std::abs(m_approximation[0]);
        while (spotToInc < m_rank) {
            bool firstPass = true;
            while(++inputSpot[spotToInc] < m_degree+1) {
                size_t spot = 0;
                for (size_t i = 0; i < m_rank; i++) {
                    spot += inputSpot[i]*multipliers[i];
                }
                m_sumAbsVal += std::abs(m_approximation[spot]);
                
                if(firstPass && spotToInc != 0) {
                    spotToInc = 0;
                }
                firstPass = false;
            }
            inputSpot[spotToInc] = 0;
            spotToInc++;
        }
        
        m_absValWasSummed = true;
    }
}


#endif /* ChebyshevApproximationND_ipp */
