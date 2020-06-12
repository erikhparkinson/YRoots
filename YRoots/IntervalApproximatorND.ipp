//
//  IntervalApproximatorND.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/9/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalApproximatorND_ipp
#define IntervalApproximatorND_ipp

template <Dimension D>
IntervalApproximator<D>::IntervalApproximator(const std::unique_ptr<FunctionInterface>& _function, size_t _approximationDegree):
m_function(_function),
m_rank(m_function->getDimension()),
m_approximationDegree(_approximationDegree),
m_sideLength(2*_approximationDegree),
m_arrayLength(power(m_sideLength, m_rank))
{
    //Alllocate memory
    m_dimensions = (int*) malloc(m_rank * sizeof (int));
    m_input = fftw_alloc_real(m_arrayLength);
    m_output= fftw_alloc_real(m_arrayLength);
    m_kinds = (fftw_r2r_kind*) malloc(m_rank * sizeof (fftw_r2r_kind));

    //Define the dimensions and kinds
    for(size_t i = 0; i < m_rank; i++) {
        m_dimensions[i] = int(m_sideLength);
        m_kinds[i] = FFTW_R2HC;
    }

    //Crete the plan
    //TODO: Decide whether to use FFTW_MEASURE or FFTW_ESTIMATE. See http://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags
    m_plan =  fftw_plan_r2r(int(m_rank), m_dimensions, m_input, m_output, m_kinds, FFTW_MEASURE | FFTW_PRESERVE_INPUT);

    //Create the precomputed points.
    m_evaluationPointsPreTransform.resize(m_arrayLength);
    m_evaluationPoints.resize(m_arrayLength);
    for(size_t i = 0; i < m_arrayLength; i++) {
        m_evaluationPointsPreTransform[i].resize(m_rank);
        m_evaluationPoints[i].resize(m_rank);
    }
    
    //Precompute the pre-transform points
    preComputeEvaluationPointsPreTransform();
    
    //Divide spots by 2
    preComputeDivideByTwoPoints();
}

template <Dimension D>
IntervalApproximator<D>::~IntervalApproximator()
{
    //Deallocate everything
    free(m_dimensions);
    fftw_free(m_input);
    fftw_free(m_output);
    free(m_kinds);
    fftw_destroy_plan(m_plan);
}

template<Dimension D>
void IntervalApproximator<D>::preComputeEvaluationPointsPreTransform()
{
    std::vector<double> chebyshevNodes;
    for(size_t i = 0; i < m_sideLength; i++) {
        chebyshevNodes.push_back(cos(i*M_PI/m_approximationDegree));
    }
    size_t divisor, number;
    for(size_t i = 0; i < m_arrayLength; i++) {
        divisor = m_arrayLength;
        number = i;
        for(size_t j = 0; j < m_rank; j++) {
            divisor /= m_sideLength;
            m_evaluationPointsPreTransform[i][j] = chebyshevNodes[number/divisor];
            number -= divisor*(number/divisor);
        }
    }
}

template<Dimension D>
void IntervalApproximator<D>::preComputeDivideByTwoPoints()
{

}


template <Dimension D>
void IntervalApproximator<D>::approximate(const Interval& _currentInterval)
{
    //Transform the evaluation points
    for(size_t j = 0; j < m_rank; j++) {
        double temp1 = _currentInterval.upperBounds[j]-_currentInterval.lowerBounds[j];
        double temp2 = _currentInterval.upperBounds[j]+_currentInterval.lowerBounds[j];
        for(size_t i = 0; i < m_evaluationPoints.size(); i++) {
            m_evaluationPoints[i][j] = (temp1*m_evaluationPointsPreTransform[i][j]+temp2)/2.0;
        }
    }
        
    //Evaluate the functions at the points
    m_function->evaluatePoints(m_evaluationPoints, m_input);
    //Divide all the inputs by degree**dimension
    size_t divisor = power(m_approximationDegree, m_rank);
    for(size_t i = 0; i < m_arrayLength; i++) {
        m_input[i] /= divisor;
    }
        
    //Take the fourier transform
    fftw_execute(m_plan);
        
    //Divide spots by 2
    for(size_t i = 0; i < m_divideByTwoPoints.size(); i++) {
        m_output[m_divideByTwoPoints[i]] /= 2;
    }
    
    printOutputArray();
}

template<Dimension D>
void IntervalApproximator<D>::printOutputArray()
{
    std::cout<<"No printing for ND!\n";
}

template<Dimension D>
void IntervalApproximator<D>::printInputArray()
{
    std::cout<<"No printing for ND!\n";
}


#endif /* IntervalApproximatorND_ipp */
