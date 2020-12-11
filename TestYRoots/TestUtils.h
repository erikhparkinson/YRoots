//
//  TestUtils.h
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 7/7/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef TestUtils_h
#define TestUtils_h

#include <iostream>

#include <cmath>
template<class T1, class T2>
bool withinEpslion(T1 a, T2 b, double epsilon = 1.e-10) {
    return std::abs(a-b) < epsilon;
}

#endif /* TestUtils_h */
