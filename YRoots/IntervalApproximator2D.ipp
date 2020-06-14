//
//  IntervalApproximator2D.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/9/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalApproximator2D_ipp
#define IntervalApproximator2D_ipp

template<>
void IntervalApproximator<Dimension::Two>::preComputeDivideByTwoPoints()
{
    for(size_t i = 0; i < m_approximationDegree; i++) {
        m_divideByTwoPoints.push_back(i);
        m_divideByTwoPoints.push_back(i + m_sideLength*m_approximationDegree);
        m_divideByTwoPoints.push_back(i*m_sideLength);
        m_divideByTwoPoints.push_back(m_approximationDegree + m_sideLength*m_approximationDegree);
    }
}

template<>
void IntervalApproximator<Dimension::Two>::printOutputArray()
{
    for(size_t i = 0; i <= m_approximationDegree; i++) {
        for(size_t j = 0; j <= m_approximationDegree; j++) {
            std::cout<<m_output[i*m_sideLength+j]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}

template<>
void IntervalApproximator<Dimension::Two>::printInputArray()
{
    for(size_t i = 0; i <= m_approximationDegree; i++) {
        for(size_t j = 0; j <= m_approximationDegree; j++) {
            std::cout<<m_input[i*m_sideLength+j]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}

#endif /* IntervalApproximator2D_ipp */
