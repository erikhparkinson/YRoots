//
//  IntervalApproximator.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/9/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalApproximator_ipp
#define IntervalApproximator_ipp

#include <fftw3.h>
#include "Function.h"
#include "utilities.h"
#include <iostream>

//TODO: Hold the function inside of here
template <Dimension D>
class IntervalApproximator
{
public:
    IntervalApproximator(const std::unique_ptr<FunctionInterface>& _function, size_t _approximationDegree);
    IntervalApproximator(IntervalApproximator const&) = delete;
    IntervalApproximator& operator=(IntervalApproximator const&) = delete;
    ~IntervalApproximator();
    
    void approximate(const Interval& _currentInterval);

    double* getOutput() {
        return m_output;
    }

    double* getInput() {
        return m_input;
    }

private:
    void preComputeEvaluationPointsPreTransform();
    void preComputeDivideByTwoPoints();
    void preComputePartialToFullTransition();
    void printOutputArray();
    void printInputArray();
    
private:
    const std::unique_ptr<FunctionInterface>& m_function;
    size_t          m_rank;
    size_t          m_approximationDegree;
    size_t          m_sideLength;
    size_t          m_arrayLength;
    int*            m_dimensions;
    double*         m_input;
    double*         m_output;
    fftw_r2r_kind*  m_kinds;
    double*         m_inputPartial;
    fftw_plan       m_plan;
    
    //For evaluating just part of the grid
    size_t          m_partialSideLength;
    size_t          m_partialArrayLength;
    std::vector<size_t> m_partialToFullTransition;
    
    //Precomputed Points
    std::vector<std::vector<double>>    m_evaluationPointsPreTransform;
    std::vector<std::vector<double>>    m_evaluationPoints;
    std::vector<size_t>                 m_divideByTwoPoints;
};

#include "IntervalApproximator1D.ipp"
#include "IntervalApproximator2D.ipp"
#include "IntervalApproximatorND.ipp"

#endif /* IntervalApproximator_ipp */
