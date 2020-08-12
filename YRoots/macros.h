//
//  macros.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 8/8/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef macros_h
#define macros_h
//Likely, Unlikely
#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)
//For Cmath
#define _USE_MATH_DEFINES

#endif /* macros_h */
