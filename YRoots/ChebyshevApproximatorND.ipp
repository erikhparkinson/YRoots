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
ChebyshevApproximator<D>::ChebyshevApproximator(size_t _rank, size_t _maxApproximationDegree, ChebyshevApproximation<D>& _approximation):
m_rank(_rank),
m_maxApproximationDegree(2*_maxApproximationDegree), //We will always double the approximation given
m_approximation(_approximation)
{
    size_t sideLength1 = m_maxApproximationDegree;
    size_t arrayLength1 = power(sideLength1, m_rank);
    size_t sideLength2 = 2*m_maxApproximationDegree;
    size_t arrayLength2 = power(sideLength2, m_rank);
    size_t partialSideLength = m_maxApproximationDegree + 1;
    size_t partialArrayLength = power(partialSideLength, m_rank);
    
    //Alllocate memory
    m_input = fftw_alloc_real(arrayLength2);
    m_output1 = fftw_alloc_real(arrayLength1);
    m_output2 = fftw_alloc_real(arrayLength2);
    m_kinds = (fftw_r2r_kind*) malloc(m_rank * sizeof (fftw_r2r_kind));
    m_inputPartial = fftw_alloc_real(partialArrayLength);

    //Define the kinds
    for(size_t i = 0; i < m_rank; i++) {
        m_kinds[i] = FFTW_R2HC;
    }
    
    //Create Interval Approximators
    size_t nextChangeDegree = _maxApproximationDegree;
    bool use1 = false;
    m_intervalApproximators.resize(m_maxApproximationDegree);
    for(size_t degree = m_maxApproximationDegree; degree > 0; degree--) {
        if(degree == nextChangeDegree) {
            nextChangeDegree /= 2;
            use1 = ~use1;
        }
        m_intervalApproximators[degree-1] = std::make_unique<IntervalApproximator<D>>(m_rank, degree, m_input, use1 ? m_output1 : m_output2, m_kinds, partialArrayLength);
    }
}

template <Dimension D>
ChebyshevApproximator<D>::~ChebyshevApproximator()
{
    //Interval approximators must be cleared first to destroy the fftw plans
    m_intervalApproximators.clear();
    
    //Deallocate everything
    fftw_free(m_inputPartial);
    free(m_kinds);
    fftw_free(m_output1);
    fftw_free(m_output2);
    fftw_free(m_input);
}


template <Dimension D>
void ChebyshevApproximator<D>::approximate(const std::unique_ptr<Function>& _function, const Interval& _currentInterval, size_t _approximationDegree)
{
    if(_approximationDegree > m_intervalApproximators.size()) {
        std::string errorMessage = "Approximation Degree is too large!";
        std::cout<<errorMessage<<"\n";
        throw std::runtime_error(errorMessage);
    }
    else if(_approximationDegree == 0) {
        std::string errorMessage = "Approximation Degree can not be 0!";
        std::cout<<errorMessage<<"\n";
        throw std::runtime_error(errorMessage);
    }
    
    m_firstApproximator = _approximationDegree-1;
    m_secondApproximator = 2*_approximationDegree-1;
    m_sideLength1 = 2*_approximationDegree;
    m_sideLength2 = 4*_approximationDegree;
    m_approxLength1 = _approximationDegree+1;
    m_approxLength2 = 2*_approximationDegree+1;

    m_intervalApproximators[m_firstApproximator]->approximate(_function, _currentInterval, false);
    m_intervalApproximators[m_secondApproximator]->approximate(_function, _currentInterval, true);
    m_infNorm = m_intervalApproximators[m_secondApproximator]->getInfoNorm();
    m_signChange = m_intervalApproximators[m_secondApproximator]->getSignChange();
    calculateApproximationError();
    
    m_approximation.setApproximation(_approximationDegree, m_sideLength1, m_intervalApproximators[m_firstApproximator]->getOutput(), m_infNorm, m_signChange, m_approximationError);
}

template <Dimension D>
void ChebyshevApproximator<D>::calculateApproximationError()
{
    double* approximation1 = m_intervalApproximators[m_firstApproximator]->getOutput();
    double* approximation2 = m_intervalApproximators[m_secondApproximator]->getOutput();
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
        while(++inputSpot[spotToInc] < m_approxLength2) {
            
            size_t spot1 = 0;
            size_t spot2 = 0;
            bool use1 = true;
            for (size_t i = 0; i < m_rank; i++) {
                use1 &= inputSpot[i] < m_approxLength1;
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

#endif /* ChebyshevApproximatorND_ipp */
