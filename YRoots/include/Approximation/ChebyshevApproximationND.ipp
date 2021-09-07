//
//  ChebyshevApproximationND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/14/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ChebyshevApproximationND_ipp
#define ChebyshevApproximationND_ipp

template <int Rank>
ChebyshevApproximation<Rank>::ChebyshevApproximation(size_t _rank):
m_rank(_rank),
m_partialSideLength(0),
m_degree(0),
m_sideLength(0),
m_approximation(nullptr),
m_infNorm(0),
m_signChange(false),
m_approximationError(0),
m_absValWasSummed(false),
m_sumAbsVal(0),
m_goodDegree(std::numeric_limits<size_t>::max() - 1)
{
    clear();
}

template <int Rank>
void ChebyshevApproximation<Rank>::setApproximation(size_t _degree, size_t _sideLength, double* _appoximation, double _infNorm, bool _signChange, double _approximationError)
{
    m_partialSideLength = _degree+1;
    m_degree = _degree*m_rank;
    m_sideLength = _sideLength;
    m_approximation = _appoximation;
    m_infNorm = _infNorm;
    m_signChange = _signChange;
    m_approximationError = _approximationError;
    clear();
    
    //Degree will never be smaller than 1 so we never need to access anything in the first two spots.
    size_t oldSize = std::max((size_t)2, m_degreeSpots.size());
    if(unlikely(m_partialSideLength + 1 > oldSize)) {
        m_degreeSpots.resize(m_partialSideLength + 1);
    }
    if(unlikely(m_degreeSpots[m_partialSideLength].size() == 0)) {
        m_degreeSpots[m_partialSideLength].resize(m_degree+1);
        setDegreeSpots(m_partialSideLength);
    }
}

template <int Rank>
void ChebyshevApproximation<Rank>::sumAbsValues() {
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
            while(++inputSpot[spotToInc] < m_partialSideLength) {
                size_t spot = 0;
                for (size_t i = 0; i < m_rank; i++) {
                    spot += inputSpot[i]*multipliers[i];
                }
                m_sumAbsVal += std::abs(m_approximation[spot]);
                
                if(spotToInc != 0) {
                    spotToInc = 0;
                }
            }
            inputSpot[spotToInc] = 0;
            spotToInc++;
        }
        
        m_absValWasSummed = true;
    }
}

template <int Rank>
inline bool ChebyshevApproximation<Rank>::isGoodApproximationSetDegree(double absApproxTol, double relApproxTol) {
    if(m_approximationError < (absApproxTol + relApproxTol*m_infNorm)) {
        m_goodDegree = m_degree;
        return true;
    }
    else{
        return false;
    }
}

template <int Rank>
inline bool ChebyshevApproximation<Rank>::isGoodApproximation(double absApproxTol, double relApproxTol) {
    return m_approximationError < (absApproxTol + relApproxTol*m_infNorm);
}

template <int Rank>
void ChebyshevApproximation<Rank>::setDegreeSpots(size_t index) {
    size_t degree = index-1;
    
    //Set up the needed variables
    std::vector<size_t> inputSpot(m_rank,0);
    std::vector<size_t> multipliers(m_rank, 1);
    for(size_t i = 1; i < m_rank; i++) {
        multipliers[i] = 2*degree*multipliers[i-1];
    }
    
    //Iterate through all the tuples of size m_rank with size up to degree.
    //Push the index of that tuple in the approximation array back into
    //m_degreeSpots at the index of the sum of the tuple.
    size_t spotToInc = 0;
    m_degreeSpots[index][0].push_back(0);
    while (spotToInc < m_rank) {
        while(++inputSpot[spotToInc] <= degree) {
            //Find the array index and sum of the tuple
            size_t arrayIndex = 0;
            size_t tupleSum = 0;
            for (size_t i = 0; i < m_rank; i++) {
                arrayIndex += inputSpot[i]*multipliers[i];
                tupleSum += inputSpot[i];
            }
            //Find the sum of the tuple
            m_degreeSpots[index][tupleSum].push_back(arrayIndex);
            
            if(spotToInc != 0) {
                spotToInc = 0;
            }
        }
        inputSpot[spotToInc] = 0;
        spotToInc++;
    }
}

template <int Rank>
bool ChebyshevApproximation<Rank>::trimCoefficients(double _absApproxTol, double _relApproxTol, size_t _targetDegree) {
    while(isGoodApproximationSetDegree(_absApproxTol, _relApproxTol)) {
        if(m_degree <= _targetDegree) {
            return true;
        }
        for(size_t& spot : m_degreeSpots[m_partialSideLength][m_degree]) {
            if(m_absValWasSummed) {
                m_sumAbsVal -= std::abs(m_approximation[spot]);
            }
            m_approximationError += std::abs(m_approximation[spot]);
            m_approximation[spot] = 0;
        }
        m_degree--;
    }
    return false;
}

template <int Rank>
void ChebyshevApproximation<Rank>::clear() {
    m_absValWasSummed = false;
    m_sumAbsVal = 0;
    m_goodDegree = std::numeric_limits<size_t>::max() - 1;
}


template <int Rank>
double* ChebyshevApproximation<Rank>::getArray() {
    return m_approximation;
}

template <int Rank>
bool ChebyshevApproximation<Rank>::isLinear() {
    return m_degree == 1;
}

template <int Rank>
double ChebyshevApproximation<Rank>::getSumAbsVal() {
    if(!m_absValWasSummed){
        sumAbsValues();
    }
    return m_sumAbsVal;
}

template <int Rank>
bool ChebyshevApproximation<Rank>::hasSignChange() {
    return m_signChange;
}

template <int Rank>
double ChebyshevApproximation<Rank>::getApproximationError() {
    return m_approximationError;
}

template <int Rank>
size_t ChebyshevApproximation<Rank>::getSideLength() {
    return m_sideLength;
}

template <int Rank>
size_t ChebyshevApproximation<Rank>::getGoodDegree() {
    return m_goodDegree;
}

template <int Rank>
size_t ChebyshevApproximation<Rank>::getDegree() {
    return m_degree;
}

template <int Rank>
size_t ChebyshevApproximation<Rank>::getPartialSideLength() {
    return m_partialSideLength;
}

#endif /* ChebyshevApproximationND_ipp */
