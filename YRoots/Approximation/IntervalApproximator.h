//
//  IntervalApproximator.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/9/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalApproximator_hpp
#define IntervalApproximator_hpp

#include <fftw3.h>
#include <iostream>
#include "../Functions/Function.h"
#include "../Utilities/utilities.h"

template <Dimension D>
class IntervalApproximator
{
public:
    IntervalApproximator(size_t _rank, size_t _approximationDegree, double* _input, double* _output, fftw_r2r_kind* _kinds, size_t _inputPartialSize);
    IntervalApproximator(IntervalApproximator const&) = delete;
    IntervalApproximator& operator=(IntervalApproximator const&) = delete;
    ~IntervalApproximator();
    
    void approximate(const Function::SharedFunctionPtr _function, const Interval& _currentInterval, bool _findInfNorm);

    bool getSignChange() {
        return m_signChange;
    }
    
    double getInfNorm() {
        return m_infNorm;
    }
    
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
    size_t          m_rank;
    size_t          m_approximationDegree;
    size_t          m_sideLength;
    size_t          m_arrayLength;
    size_t          m_partialSideLength;
    size_t          m_partialArrayLength;
    
    int*            m_dimensions;
    double*         m_input;
    double*         m_output;
    fftw_r2r_kind*  m_kinds;
    std::vector<double> m_inputPartial;
    fftw_plan       m_plan;
    
    //For evaluating just part of the grid
    std::vector<size_t> m_partialToFullTransition;
    
    //Precomputed Points
    std::vector<std::vector<double> >    m_evaluationPointsPreTransform;
    std::vector<std::vector<double> >    m_evaluationPoints;
    std::vector<size_t>                 m_divideByTwoPoints;
    
    //Other
    double          m_infNorm;
    bool            m_signChange;
    
    static size_t           m_timerInitIndex;
    static size_t           m_timerPlanIndex;
    static size_t           m_timerIntervalApproximatorIndex;
    static size_t           m_timerFFT;
    Timer&                  m_timer = Timer::getInstance();
};

template<Dimension D>
size_t IntervalApproximator<D>::m_timerIntervalApproximatorIndex = -1;
template<Dimension D>
size_t IntervalApproximator<D>::m_timerFFT = -1;
template<Dimension D>
size_t IntervalApproximator<D>::m_timerInitIndex = -1;
template<Dimension D>
size_t IntervalApproximator<D>::m_timerPlanIndex = -1;

#include "IntervalApproximator1D.ipp"
#include "IntervalApproximator2D.ipp"
#include "IntervalApproximatorND.ipp"

#endif /* IntervalApproximator_hpp */
