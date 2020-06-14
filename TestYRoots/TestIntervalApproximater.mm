//
//  TestIntervalApproximater.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 6/13/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef TestIntervalApproximater_mm
#define TestIntervalApproximater_mm

#import <XCTest/XCTest.h>
#include "IntervalApproximator.h"
#include "PowerBasisPolynomial.h"

@interface TestIntervalApproximater : XCTestCase

@end

@implementation TestIntervalApproximater

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

template<class T>
bool withinEpslion(T a, T b, double epsilon = 1.e-10) {
    return std::abs(a-b) < epsilon;
}

- (void)testBasic1D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    std::string functionString = "5+x1^2";
    size_t approximationDegree = 3;
    
    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    IntervalApproximator<Dimension::One> intervalApproximator(function, approximationDegree);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0);
    intervalApproximator.approximate(currentInterval);

    XCTAssert(withinEpslion(intervalApproximator.getOutput()[0], 5.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[1], 0.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[2], 0.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[3], 0.0));
}

- (void)testBasic2D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    variablesNames.push_back("x2");
    std::string functionString = "5+x1^2+x2";
    size_t approximationDegree = 3;
    
    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    IntervalApproximator<Dimension::Two> intervalApproximator(function, approximationDegree);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0);
    intervalApproximator.approximate(currentInterval);

    XCTAssert(withinEpslion(intervalApproximator.getOutput()[0], 5.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[1], 0.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[2], 0.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[3], 0.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[6], 1.0));
    for(size_t i = 7;  i < 10; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
    for(size_t i = 12;  i < 16; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
    for(size_t i = 18;  i < 22; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
}

- (void)testBasic3D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    variablesNames.push_back("x2");
    variablesNames.push_back("x3");
    std::string functionString = "5+x1^2+x2+5*x1*x2*x3";
    size_t approximationDegree = 2;
    
    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    IntervalApproximator<Dimension::Three> intervalApproximator(function, approximationDegree);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0);
    intervalApproximator.approximate(currentInterval);
    
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[0], 5.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[1], 0.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[2], 0.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[4], 1.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[5], 0.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[6], 0.0));
    for(size_t i = 8;  i < 10; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
    
    for(size_t i = 16;  i < 18; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[20], 0.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[21], 5.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[22], 0.0));
    for(size_t i = 24;  i < 26; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
    
    for(size_t i = 32;  i < 36; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
    for(size_t i = 36;  i < 38; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
    for(size_t i = 40;  i < 42; i++) {
        XCTAssert(withinEpslion(intervalApproximator.getOutput()[i], 0.0));
    }
}

@end

#endif
