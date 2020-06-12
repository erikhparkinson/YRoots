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
ChebyshevApproximator<D>::ChebyshevApproximator(const std::unique_ptr<FunctionInterface>& _function, size_t approximationDegree):
m_intervalApproximator1(_function, approximationDegree),
m_intervalApproximator2(_function, 2*approximationDegree)
{
    
}

template <Dimension D>
ChebyshevApproximator<D>::~ChebyshevApproximator()
{

}


template <Dimension D>
void ChebyshevApproximator<D>::approximate(const Interval& _currentInterval)
{
    m_intervalApproximator1.approximate(_currentInterval);
    m_intervalApproximator2.approximate(_currentInterval);
}



#endif /* ChebyshevApproximatorND_ipp */
