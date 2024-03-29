//
//  IntervalApproximatorND.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/9/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalApproximatorND_ipp
#define IntervalApproximatorND_ipp

template <int Rank>
IntervalApproximator<Rank>::IntervalApproximator(size_t _rank, size_t _approximationDegree, double* _input, double* _output, fftw_r2r_kind* _kinds, size_t _inputPartialSize):
m_rank(_rank),
m_approximationDegree(_approximationDegree),
m_sideLength(2*_approximationDegree),
m_arrayLength(power(m_sideLength, m_rank)),
m_partialSideLength(_approximationDegree + 1),
m_partialArrayLength(power(m_partialSideLength, m_rank)),
m_input(_input),
m_output(_output),
m_kinds(_kinds),
m_infNorm(0),
m_signChange(false)
{
    m_inputPartial.resize(_inputPartialSize);
    
    //Alllocate memory
    m_dimensions = (int*) malloc(m_rank * sizeof (int));
    
    //Define the dimensions
    for(size_t i = 0; i < m_rank; i++) {
        m_dimensions[i] = int(m_sideLength);
    }

    //Crete the plan
    //Options FFTW_EXHAUSTIVE, FFTW_PATIENT, FFTW_MEASURE, FFTW_ESTIMATE. See http://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags.
    m_plan =  fftw_plan_r2r(int(m_rank), m_dimensions, m_input, m_output, m_kinds, FFTW_PATIENT | FFTW_DESTROY_INPUT);

    //Create the precomputed points.
    m_evaluationPointsPreTransform.resize(m_rank);
    m_evaluationPoints.resize(m_rank);
    for(size_t i = 0; i < m_rank; i++) {
        m_evaluationPointsPreTransform[i].resize(m_partialSideLength);
        m_evaluationPoints[i].resize(m_partialSideLength);
    }
    
    //Precompute the pre-transform points
    preComputeEvaluationPointsPreTransform();
    
    //Divide spots by 2
    preComputeDivideByTwoPoints();
    
    //Precompute the partial to full transition
    preComputePartialToFullTransition();
        
    m_timer.registerTimer(m_timerIntervalApproximatorIndex, "Interval Approximator");
    m_timer.registerTimer(m_timerFFT, "FFT");
    m_timer.registerTimer(m_timerEvalGrid, "Evaluate Grid");
}

template <int Rank>
IntervalApproximator<Rank>::~IntervalApproximator()
{
    //Deallocate everything
    fftw_destroy_plan(m_plan);
    free(m_dimensions);
}

template<int Rank>
void IntervalApproximator<Rank>::preComputeEvaluationPointsPreTransform()
{
    for(size_t i = 0; i < m_partialSideLength; i++) {
        double val = cos(i*M_PI/m_approximationDegree);
        for(size_t j = 0; j < m_rank; j++) {
            m_evaluationPointsPreTransform[j][i] = val;
        }
    }
}

template<int Rank>
void IntervalApproximator<Rank>::preComputeDivideByTwoPoints()
{
    size_t divisor, number;
    size_t fixedNums = power(m_approximationDegree+1, m_rank-1);
    
    std::vector<size_t> powers;
    size_t temp = 1;
    for(size_t i = 0; i < m_rank; i++) {
        powers.push_back(temp);
        temp *= m_sideLength;
    }
    
    std::vector<size_t> spots;
    spots.resize(m_rank);

    for(size_t fixedVar = 0; fixedVar < m_rank; fixedVar++) {
        for(size_t i = 0; i < fixedNums; i++) {
            divisor = fixedNums;
            number = i;
            for(size_t j = 0; j < m_rank; j++) {
                if(j != fixedVar) {
                    divisor /= (m_approximationDegree+1);
                    spots[j] = number/divisor;
                    number -= divisor*(number/divisor);
                }
            }
            
            //Add it
            size_t result = 0;
            spots[fixedVar] = 0;
            for(size_t j =0; j < m_rank; j++) {
                result += spots[j] * powers[j];
            }
            m_divideByTwoPoints.push_back(result);
            
            result = 0;
            spots[fixedVar] = m_approximationDegree;
            for(size_t j =0; j < m_rank; j++) {
                result += spots[j] * powers[j];
            }
            m_divideByTwoPoints.push_back(result);
        }
    }
}

template<int Rank>
void IntervalApproximator<Rank>::preComputePartialToFullTransition()
{
    //Set up the needed variables
    std::vector<size_t> inputSpot(m_rank);
    
    //Iterate through all the combinations
    size_t spotToInc = 0;
    m_partialToFullTransition.push_back(0);
    while (spotToInc < m_rank) {
        while(++inputSpot[spotToInc] < m_sideLength) {
            size_t temp = 1;
            size_t result = 0;
            for (size_t i = 0; i < m_rank; i++) {
                result += (inputSpot[i] <= m_approximationDegree ? inputSpot[i] : m_sideLength - inputSpot[i])*temp;
                temp *= (m_approximationDegree+1);
            }
            m_partialToFullTransition.push_back(result);
            if(spotToInc != 0) {
                spotToInc = 0;
            }
        }
        inputSpot[spotToInc] = 0;
        spotToInc++;
    }
}

template <int Rank>
void IntervalApproximator<Rank>::approximate(const Function::SharedFunctionPtr _function, const Interval& _currentInterval, bool _findInfNorm)
{
    m_timer.startTimer(m_timerIntervalApproximatorIndex);

    //Transform the evaluation points
    for(size_t j = 0; j < m_rank; j++) {
        double temp1 = _currentInterval.upperBounds[j] - _currentInterval.lowerBounds[j];
        double temp2 = _currentInterval.upperBounds[j] + _currentInterval.lowerBounds[j];
        for(size_t i = 0; i < m_partialSideLength; i++) {
            m_evaluationPoints[j][i] = (temp1*m_evaluationPointsPreTransform[j][i] + temp2)/2.0;
        }
    }
            
    //Evaluate the functions at the points
    double divisor = static_cast<double>(power(m_approximationDegree, m_rank));
    m_timer.startTimer(m_timerEvalGrid);
    _function->evaluateGrid(m_evaluationPoints, m_inputPartial);
    m_timer.stopTimer(m_timerEvalGrid);

    if(_findInfNorm) {
        m_infNorm = 0;
        bool hasPos = false;
        bool hasNeg = false;
        for(size_t i = 0; i < m_partialArrayLength; i++) {
            m_infNorm = std::max(m_infNorm, std::abs(m_inputPartial[i]));
            hasPos |= m_inputPartial[i] > 0;
            hasNeg |= m_inputPartial[i] < 0;
        }
        m_signChange = hasPos && hasNeg;
    }
    
    //Fill in the full input
    for(size_t i = 0; i < m_arrayLength; i++) {
        m_input[i] = m_inputPartial[m_partialToFullTransition[i]] / divisor;
    }
    
    //Take the fourier transform
    m_timer.startTimer(m_timerFFT);
    fftw_execute(m_plan);
    m_timer.stopTimer(m_timerFFT);

    //Divide spots by 2
    for(size_t i = 0; i < m_divideByTwoPoints.size(); i++) {
        m_output[m_divideByTwoPoints[i]] /= 2;
    }
    
    //Print input and output. TODO: Remove when I'm confident things work.
//    std::cout<<"Input:\n";
//    printInputArray();
//    std::cout<<"Output:\n";
//    printOutputArray();
    
    m_timer.stopTimer(m_timerIntervalApproximatorIndex);
}

template<int Rank>
void IntervalApproximator<Rank>::printOutputArray()
{
    std::cout<<"No printing for ND!\n";
}

template<int Rank>
void IntervalApproximator<Rank>::printInputArray()
{
    std::cout<<"No printing for ND!\n";
}


#endif /* IntervalApproximatorND_ipp */
