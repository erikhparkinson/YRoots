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
    m_evaluationPointsPreTransform.resize(m_rank);
    m_evaluationPoints.resize(m_rank);
    for(size_t i = 0; i < m_rank; i++) {
        m_evaluationPointsPreTransform[i].resize(m_sideLength);
        m_evaluationPoints[i].resize(m_sideLength);
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
    for(size_t i = 0; i < m_sideLength; i++) {
        double val = cos(i*M_PI/m_approximationDegree);
        for(size_t j = 0; j < m_rank; j++) {
            m_evaluationPointsPreTransform[j][i] = val;
        }
    }
}

template<Dimension D>
void IntervalApproximator<D>::preComputeDivideByTwoPoints()
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


template <Dimension D>
void IntervalApproximator<D>::approximate(const Interval& _currentInterval)
{
    //Transform the evaluation points
    for(size_t j = 0; j < m_rank; j++) {
        double temp1 = _currentInterval.upperBounds[j]-_currentInterval.lowerBounds[j];
        double temp2 = _currentInterval.upperBounds[j]+_currentInterval.lowerBounds[j];
        for(size_t i = 0; i < m_sideLength; i++) {
            m_evaluationPoints[j][i] = (temp1*m_evaluationPointsPreTransform[j][i]+temp2)/2.0;
        }
    }
            
    //Evaluate the functions at the points
    m_function->evaluateGrid(m_evaluationPoints, m_input);
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
