//
//  ChebyshevApproximator.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/5/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ChebyshevApproximator_h
#define ChebyshevApproximator_h

#include "Approximation/IntervalApproximator.hpp"
#include "Approximation/ChebyshevApproximation.hpp"
#include "Utilities/Timer.hpp"

template <int Rank>
class ChebyshevApproximator
{
public:
    ChebyshevApproximator(size_t _rank, size_t _maxApproximationDegree, ChebyshevApproximation<Rank>& _approximation);
    ~ChebyshevApproximator();
    
    void approximate(const Function::SharedFunctionPtr _function, const Interval& _currentInterval, size_t _approximationDegree);
    double getAbsApproxTol(const Function::SharedFunctionPtr _function, const Interval& _currentInterval);

    bool hasSignChange() {
        return m_signChange;
    }
    
    ChebyshevApproximation<Rank>& getApproximation() {
        return m_approximation;
    }
    
private:
    void calculateApproximationError();
    
private:
    size_t                                  m_rank;
    size_t                                  m_maxApproximationDegree;
    
    double*                                 m_input;
    double*                                 m_output1;
    double*                                 m_output2;
    fftw_r2r_kind*                          m_kinds;
    double*                                 m_inputPartial;
    
    std::vector<std::unique_ptr<IntervalApproximator<Rank>>>    m_intervalApproximators;
    size_t                                  m_firstApproximator;
    size_t                                  m_secondApproximator;
    size_t                                  m_sideLength1;
    size_t                                  m_sideLength2;
    size_t                                  m_approxLength1;
    size_t                                  m_approxLength2;

    double                                  m_infNorm;
    bool                                    m_signChange;
    double                                  m_approximationError;
    
    Interval                                m_absApproxErrorCalcInterval;
    
    ChebyshevApproximation<Rank>&              m_approximation;
    
    static size_t           m_timerFullApproximateIndex;
    static size_t           m_timerAbsApproxErrorCalcIndex;
    Timer&                  m_timer = Timer::getInstance();
};

template<int Rank>
size_t ChebyshevApproximator<Rank>::m_timerFullApproximateIndex = -1;
template<int Rank>
size_t ChebyshevApproximator<Rank>::m_timerAbsApproxErrorCalcIndex = -1;

    
#include "ChebyshevApproximatorND.ipp"

#endif /* ChebyshevApproximator_h */
