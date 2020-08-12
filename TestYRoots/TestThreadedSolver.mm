//
//  TestThreadedSolver.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 7/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "TestUtils.h"
#include "ThreadedSolver.h"

@interface TestThreadedSolver : XCTestCase

@end

@implementation TestThreadedSolver

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testBasic1D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    std::string functionString = "-1+2*x1^25";
    
    std::unique_ptr<FunctionInterface> function = std::make_unique<PowerBasisPolynomial>(functionString, variablesNames);
    Interval startInterval;
    startInterval.lowerBounds.push_back(-1.0);
    startInterval.upperBounds.push_back(1.0);
    
    for(size_t numThreads = 1; numThreads <= 4; numThreads++) {
        
    std::vector<std::vector<std::unique_ptr<FunctionInterface>>> allFunctions;
    allFunctions.resize(numThreads);
    for(size_t i = 0; i < numThreads; i++) {
        allFunctions[i].emplace_back(std::make_unique<PowerBasisPolynomial>(functionString, variablesNames));
    }
    ThreadedSolver<Dimension::One> solver(allFunctions, numThreads, startInterval);
        
    size_t trials = 1;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        solver.solve();
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "Solve with " << numThreads << " threads takes " <<nanos/(trials * 1000000)<< "ms.\n";

    XCTAssert(true);
    }
}

@end
