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
    //They all use the same input vector, the output is the same
    size_t nextChangeDegree = _maxApproximationDegree;
    bool use1 = false;
    m_intervalApproximators.resize(m_maxApproximationDegree);
    for(size_t degree = m_maxApproximationDegree; degree > 0; degree--) {
        if(degree == nextChangeDegree) {
            nextChangeDegree /= 2;
            use1 = !use1;
        }
        m_intervalApproximators[degree-1] = std::make_unique<IntervalApproximator<D>>(m_rank, degree, m_input, use1 ? m_output1 : m_output2, m_kinds, partialArrayLength);
    }
    
    //Initialize m_absApproxErrorCalcInterval
    m_absApproxErrorCalcInterval.lowerBounds.resize(m_rank);
    m_absApproxErrorCalcInterval.upperBounds.resize(m_rank);
    
    m_timer.registerTimer(m_timerFullApproximateIndex, "Cheb Approximator Full");
    m_timer.registerTimer(m_timerAbsApproxErrorCalcIndex, "Cheb Approximator Abs Approx Estimate");
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
void ChebyshevApproximator<D>::approximate(const Function::SharedFunctionPtr _function, const Interval& _currentInterval, size_t _approximationDegree)
{
    m_timer.startTimer(m_timerFullApproximateIndex);
    if(_approximationDegree > m_intervalApproximators.size()) {
        printAndThrowRuntimeError("Approximation Degree is too large!");
    }
    else if(_approximationDegree == 0) {
        printAndThrowRuntimeError("Approximation Degree can not be 0!");
    }
    
    m_firstApproximator = _approximationDegree-1;
    m_secondApproximator = 2*_approximationDegree-1;
    m_sideLength1 = 2*_approximationDegree;
    m_sideLength2 = 4*_approximationDegree;
    m_approxLength1 = _approximationDegree+1;
    m_approxLength2 = 2*_approximationDegree+1;

    m_intervalApproximators[m_firstApproximator]->approximate(_function, _currentInterval, false);
    m_intervalApproximators[m_secondApproximator]->approximate(_function, _currentInterval, true);
    m_infNorm = m_intervalApproximators[m_secondApproximator]->getInfNorm();
    m_signChange = m_intervalApproximators[m_secondApproximator]->getSignChange();
    calculateApproximationError();
        
    m_approximation.setApproximation(_approximationDegree, m_sideLength1, m_intervalApproximators[m_firstApproximator]->getOutput(), m_infNorm, m_signChange, m_approximationError);
    
    m_timer.stopTimer(m_timerFullApproximateIndex);
}

template <Dimension D>
double ChebyshevApproximator<D>::getAbsApproxTol(const Function::SharedFunctionPtr _function, const Interval& _currentInterval, size_t _approximationDegree)
{
    m_timer.startTimer(m_timerAbsApproxErrorCalcIndex);

    m_firstApproximator = _approximationDegree-1;
    m_secondApproximator = 2*_approximationDegree-1;
    m_sideLength1 = 2*_approximationDegree;
    m_sideLength2 = 4*_approximationDegree;
    m_approxLength1 = _approximationDegree+1;
    m_approxLength2 = 2*_approximationDegree+1;
    
    //Get an interval
    const double linearization_size = 1e-14; //Make sure this isn't too small. Otherwise, we could get an interval where everything in it evaluates to the same point.
    for(size_t i = 0; i < m_rank; i++) {
        const double rand = 0.51234127384517283654; //TODO: Rand uniform [0,1]. Make sure the random point is not 0!
        const double randomPoint = _currentInterval.lowerBounds[i] + rand * (_currentInterval.upperBounds[i] - _currentInterval.lowerBounds[i]);
        m_absApproxErrorCalcInterval.lowerBounds[i] = randomPoint * (1+linearization_size);
        m_absApproxErrorCalcInterval.upperBounds[i] = randomPoint * (1-linearization_size);
    }

    m_intervalApproximators[m_firstApproximator]->approximate(_function, m_absApproxErrorCalcInterval, false);
    m_intervalApproximators[m_secondApproximator]->approximate(_function, m_absApproxErrorCalcInterval, false);
    calculateApproximationError();
    
    m_timer.stopTimer(m_timerAbsApproxErrorCalcIndex);

    return m_approximationError*10;
}

template <Dimension D>
void ChebyshevApproximator<D>::calculateApproximationError()
{
    double* approximation1 = m_intervalApproximators[m_firstApproximator]->getOutput();
    double* approximation2 = m_intervalApproximators[m_secondApproximator]->getOutput();
    assert(approximation1 != approximation2);
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
                m_approximationError += std::abs(approximation2[spot2] - approximation1[spot1]);
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
