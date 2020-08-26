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
    ChebyshevApproximation<Dimension::One> chebApproximation(1);
    ChebyshevApproximator<Dimension::One> chebyshevApproximator(1, approximationDegree, chebApproximation);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0);
    chebyshevApproximator.approximate(function, currentInterval, approximationDegree);

    double* approximation = chebApproximation.getArray();
    
    XCTAssert(withinEpslion(approximation[0], 5.5));
    XCTAssert(withinEpslion(approximation[1], 0.0));
    XCTAssert(withinEpslion(approximation[2], 0.5));
    XCTAssert(withinEpslion(approximation[3], 0.0));
    
    XCTAssertTrue(chebApproximation.isGoodApproximation(1e-15, 1e-15));
}

- (void)testBasic2D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    variablesNames.push_back("x2");
    std::string functionString = "5+x1^2+x2";
    size_t approximationDegree = 3;
    
    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    ChebyshevApproximation<Dimension::Two> chebApproximation(2);
    ChebyshevApproximator<Dimension::Two> chebyshevApproximator(2, approximationDegree, chebApproximation);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0);
    chebyshevApproximator.approximate(function, currentInterval, approximationDegree);

    double* approximation = chebApproximation.getArray();
        
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
    
    XCTAssertTrue(chebApproximation.isGoodApproximation(1e-15, 1e-15));
}

- (void)testTimingTemp {
    size_t rank = 2;
    size_t approximationDegree = 5;    
    size_t degreePoly = 40;

    std::vector<std::string> variablesNames;
    std::string functionString = "1+";
    for(size_t i = 0; i < rank; i++) {
        variablesNames.push_back("x" + std::to_string(i));
        functionString += "x" + std::to_string(i) + "^" + std::to_string(degreePoly);
        if (i+1 != rank) {
            functionString += "*";
        }
    }

    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    ChebyshevApproximation<Dimension::Two> chebApproximation(2);
    ChebyshevApproximator<Dimension::Two> chebyshevApproximator(rank, approximationDegree, chebApproximation);
    Interval currentInterval;
    for(size_t i = 0; i < rank; i++) {
        currentInterval.lowerBounds.push_back(-1.0);
        currentInterval.upperBounds.push_back(1.0);
    }

    size_t trials = 1000;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        chebyshevApproximator.approximate(function, currentInterval, approximationDegree);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "\nCreating a " << approximationDegree << " degree chebyshev approximaiton for a Degree " << degreePoly << " Dimension " << rank << " Polynomial takes " <<nanos/(trials*1000)<< "us.\n\n";
}


@end
