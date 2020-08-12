//
//  TestMain.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 7/7/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "TestPolynomialPowerBasis.mm"
#include "TestIntervalApproximater.mm"
#include "TestChebyshevApproximator.mm"
#include "TestThreadedSolver.mm"
#include "TestConcurrent.mm"
#include "TestChebyshevApproximation.mm"
#include "TestIntervalChecker.mm"

@interface TestMain : XCTestCase

@end

@implementation TestMain

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

@end
