//
//  ChebyshevApproximatorND.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/8/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ChebyshevApproximatorND_ipp
#define ChebyshevApproximatorND_ipp

template <Dimension D>
ChebyshevApproximator<D>::ChebyshevApproximator(const std::unique_ptr<FunctionInterface>& _function, size_t _approximationDegree):
m_intervalApproximator1(_function, _approximationDegree),
m_intervalApproximator2(_function, 2*_approximationDegree),
m_rank(_function->getDimension()),
m_approximationDegree(_approximationDegree),
m_sideLength1(m_approximationDegree*2),
m_sideLength2(m_sideLength1*2)
{
    
}

template <Dimension D>
ChebyshevApproximator<D>::~ChebyshevApproximator()
{

}


template <Dimension D>
void ChebyshevApproximator<D>::approximate(const Interval& _currentInterval)
{
    m_intervalApproximator1.approximate(_currentInterval, false);
    m_intervalApproximator2.approximate(_currentInterval, true);
    m_infNorm = m_intervalApproximator2.getInfoNorm();
    m_signChange = m_intervalApproximator2.getSignChange();
    calculateApproximationError();
}

template <Dimension D>
void ChebyshevApproximator<D>::calculateApproximationError()
{
    double* approximation1 = m_intervalApproximator1.getOutput();
    double* approximation2 = m_intervalApproximator2.getOutput();
    m_approximationError = 0.0;
    
    //Set up the needed variables
    std::vector<size_t> inputSpot(m_rank,0);
    std::vector<size_t> multipliers1(m_rank, 1);
    std::vector<size_t> multipliers2(m_rank, 1);
    for(size_t i = 1; i < m_rank; i++) {
        multipliers1[i] = m_sideLength1*multipliers1[i-1];
        multipliers2[i] = m_sideLength2*multipliers2[i-1];
    }
    
    //Iterate through all the combinations
    size_t spotToInc = 0;
    m_approximationError += std::abs(approximation1[0] - approximation2[0]);
    while (spotToInc < m_rank) {
        bool firstPass = true;
        while(++inputSpot[spotToInc] < m_sideLength2) {
            
            size_t spot1 = 0;
            size_t spot2 = 0;
            bool use1 = true;
            for (size_t i = 0; i < m_rank; i++) {
                use1 &= inputSpot[i] < m_sideLength1;
                spot1 += inputSpot[i]*multipliers1[i];
                spot2 += inputSpot[i]*multipliers2[i];
            }
            if(use1) {
                m_approximationError += std::abs(approximation2[spot2] - (approximation1[spot1]));
            }
            else {
                m_approximationError += std::abs(approximation2[spot2]);
            }
            
            if(firstPass && spotToInc != 0) {
                spotToInc = 0;
            }
            firstPass = false;
        }
        inputSpot[spotToInc] = 0;
        spotToInc++;
    }
}


template <Dimension D>
bool ChebyshevApproximator<D>::isGoodApproximation(SubdivisionParameters& _subdivisionParameters) {
    return m_approximationError < (_subdivisionParameters.absApproxTol + _subdivisionParameters.relApproxTol*m_infNorm);
}


#endif /* ChebyshevApproximatorND_ipp */
