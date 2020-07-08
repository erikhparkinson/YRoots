//
//  TestChebyshevApproximator.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 6/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "TestUtils.h"
#include "ChebyshevApproximator.h"

@interface TestChebyshevApproximator : XCTestCase

@end

@implementation TestChebyshevApproximator

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testBasic1D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    std::string functionString = "5+x1^2";
    size_t approximationDegree = 3;
    
    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    ChebyshevApproximator<Dimension::One> chebyshevApproximator(function, approximationDegree);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0);
    chebyshevApproximator.approximate(currentInterval);

    double* approximation = chebyshevApproximator.getApproximation();
    
    XCTAssert(withinEpslion(approximation[0], 5.5));
    XCTAssert(withinEpslion(approximation[1], 0.0));
    XCTAssert(withinEpslion(approximation[2], 0.5));
    XCTAssert(withinEpslion(approximation[3], 0.0));
    
    SubdivisionParameters subdivisionParameters;
    subdivisionParameters.relApproxTol = 1e-15;
    subdivisionParameters.absApproxTol = 1e-15;
    XCTAssertTrue(chebyshevApproximator.isGoodApproximation(subdivisionParameters));
}

- (void)testBasic2D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    variablesNames.push_back("x2");
    std::string functionString = "5+x1^2+x2";
    size_t approximationDegree = 3;
    
    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    ChebyshevApproximator<Dimension::Two> chebyshevApproximator(function, approximationDegree);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0);
    chebyshevApproximator.approximate(currentInterval);

    double* approximation = chebyshevApproximator.getApproximation();
    
    XCTAssert(withinEpslion(approximation[0], 5.5));
    XCTAssert(withinEpslion(approximation[1], 0.0));
    XCTAssert(withinEpslion(approximation[2], 0.5));
    XCTAssert(withinEpslion(approximation[3], 0.0));
    XCTAssert(withinEpslion(approximation[6], 1.0));
    for(size_t i = 7;  i < 10; i++) {
        XCTAssert(withinEpslion(approximation[i], 0.0));
    }
    for(size_t i = 12;  i < 16; i++) {
        XCTAssert(withinEpslion(approximation[i], 0.0));
    }
    for(size_t i = 18;  i < 22; i++) {
        XCTAssert(withinEpslion(approximation[i], 0.0));
    }
    
    SubdivisionParameters subdivisionParameters;
    subdivisionParameters.relApproxTol = 1e-15;
    subdivisionParameters.absApproxTol = 1e-15;
    XCTAssertTrue(chebyshevApproximator.isGoodApproximation(subdivisionParameters));
}

@end
