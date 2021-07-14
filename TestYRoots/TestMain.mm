//
//  TestMain.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 7/7/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef USE_TIMING
#define USE_TIMING
#endif

#define TESTING

#import <XCTest/XCTest.h>
#include "FFTW/include/fftw3.h"

@interface TestMain : XCTestCase

@end

@implementation TestMain

//TODO: Find a better way to do this
//Variables for the tests.
size_t              m_rank;
size_t              m_approximationDegree;
size_t              m_sideLength;
size_t              m_arrayLength;
size_t              m_partialSideLength;
size_t              m_partialArrayLength;

bool                m_allocated;
double*             m_input;
double*             m_output;
fftw_r2r_kind*      m_kinds;
double*             m_inputPartial;
double*             m_approximation;

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

@end

#include "TestPolynomialPowerBasis.mm"
#include "TestIntervalApproximater.mm"
#include "TestChebyshevApproximator.mm"
#include "TestThreadedSolver.mm"
#include "TestConcurrent.mm"
#include "TestChebyshevApproximation.mm"
#include "TestIntervalChecker.mm"
#include "TestFunction.mm"
#include "TestQuadraticCheck.mm"
