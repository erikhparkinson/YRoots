//
//  IntervalApproximator1D.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/9/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalApproximator1D_ipp
#define IntervalApproximator1D_ipp

template<>
void IntervalApproximator<1>::preComputeDivideByTwoPoints()
{
    m_divideByTwoPoints.push_back(0);
    m_divideByTwoPoints.push_back(m_approximationDegree);
}

template<>
void IntervalApproximator<1>::printOutputArray()
{
    for(size_t i = 0; i <= m_approximationDegree; i++) {
        std::cout<<m_output[i]<<" ";
    }
    std::cout<<"\n";
}

template<>
void IntervalApproximator<1>::printInputArray()
{
    for(size_t i = 0; i <= m_approximationDegree; i++) {
        std::cout<<m_input[i]<<" ";
    }
    std::cout<<"\n";
}


#endif /* IntervalApproximator1D_ipp */
