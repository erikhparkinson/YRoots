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
    //Disable the timer as default, as running with multiple threads is an issue
    Timer::disable();
}

- (void)tearDown {
}

- (void)testBasic1D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    std::vector<std::string> subfunctionNames;
    std::string functionString = "-1+2*x1^20";
    
    Interval startInterval;
    startInterval.lowerBounds.push_back(-1.0);
    startInterval.upperBounds.push_back(1.0);
    
    std::cout<<"\n";
    for(size_t numThreads = 1; numThreads <= 4; numThreads++) {
        std::vector<std::unique_ptr<Function>> functions;
        functions.emplace_back(std::make_unique<Function>(functionString, variablesNames, subfunctionNames));
        ThreadedSolver<Dimension::One> solver(functions, numThreads, startInterval);
            
        size_t trials = 1;
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        for(size_t i = 0; i < trials; i++) {
            solver.solve();
        }
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

        double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
        std::cout << "Solve with " << numThreads << " threads takes " <<nanos/(trials * 1000)<< "us.\n";

        std::vector<FoundRoot> foundRoots =  solver.getRoots();
        for(FoundRoot& root : foundRoots) {
            XCTAssert(withinEpslion(power(root.root[0], 20), 0.5));
        }
    }
    std::cout<<"\n";
}

- (void)testHarder1D {
    //Parameters
    std::vector<double> upperBounds;
    std::vector<uint16_t> powerNums;

    //Set up the tests
    upperBounds.push_back(10000.0); powerNums.push_back(1);
    upperBounds.push_back(100.0); powerNums.push_back(2);
    upperBounds.push_back(21.0); powerNums.push_back(3);
    upperBounds.push_back(10.0); powerNums.push_back(4);

    std::cout<<"\n";
    for(size_t testNum = 0; testNum < upperBounds.size(); testNum++) {
        //Get the variables
        double upperBound = upperBounds[testNum];
        uint16_t powerNum = powerNums[testNum];
        size_t numThreads = 1;

        std::vector<std::string> variablesNames;
        variablesNames.push_back("x");
        std::vector<std::string> subfunctionNames;
        //sin(x^powerNum)
        //Roots at x^powerNum = k*pi
        //There are upperBound^powerNum/pi roots
        std::string functionString = "sin(x^" + std::to_string(powerNum) + ")";
            
        //Declare the intervals
        Interval startInterval;
        startInterval.lowerBounds.push_back(0.5);
        startInterval.upperBounds.push_back(upperBound);

        //Set up the solver
        std::vector<std::unique_ptr<Function>> functions;
        functions.emplace_back(std::make_unique<Function>(functionString, variablesNames, subfunctionNames));
        ThreadedSolver<Dimension::One> solver(functions, numThreads, startInterval);
            
        //Solve it
        size_t trials = 1;
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        for(size_t i = 0; i < trials; i++) {
            solver.solve();
        }
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
        double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
        std::cout << "Solve " << functionString << " on "<< startInterval.toString() << " with " << numThreads << " threads takes " <<nanos/(trials * 1e6)<< "ms.\n";

        //Assert the solutions
        std::vector<FoundRoot> foundRoots =  solver.getRoots();
        size_t expectedRoots = power(upperBound, powerNum) / M_PI;
        for(FoundRoot& root : foundRoots) {
            double rootModPi = fmod(power(root.root[0], powerNum), M_PI);
            XCTAssert(withinEpslion(rootModPi, 0.0) || withinEpslion(rootModPi, M_PI));
        }
        XCTAssert(expectedRoots == foundRoots.size());
    }
    std::cout<<"\n";
}

- (void)testHarder2D {
    Timer::enable();

    //Parameters
    std::vector<double> upperBoundsX;
    std::vector<double> upperBoundsY;
    std::vector<uint16_t> powerNumsX;
    std::vector<uint16_t> powerNumsY;
    std::vector<double> powerNumsX2;
    std::vector<double> powerNumsY2;

    //First Test
    upperBoundsX.push_back(50); powerNumsX.push_back(1); powerNumsX2.push_back(1.0);
    upperBoundsY.push_back(50); powerNumsY.push_back(1); powerNumsY2.push_back(1.0);
    //Second Test
    upperBoundsX.push_back(50); powerNumsX.push_back(1); powerNumsX2.push_back(1.2);
    upperBoundsY.push_back(50); powerNumsY.push_back(1); powerNumsY2.push_back(0.8);
    //Third Test
    upperBoundsX.push_back(50); powerNumsX.push_back(1.5); powerNumsX2.push_back(1.2);
    upperBoundsY.push_back(50); powerNumsY.push_back(1.3); powerNumsY2.push_back(0.8);

    std::vector<std::string> variablesNames;
    variablesNames.push_back("x");
    variablesNames.push_back("y");
    std::vector<std::string> subfunctionNames;
    
    std::cout<<"\n";
    for(size_t testNum = 0; testNum < upperBoundsX.size(); testNum++) {
        //Get the variables
        double upperBoundX = upperBoundsX[testNum];
        double upperBoundY = upperBoundsY[testNum];
        uint16_t powerNumX = powerNumsX[testNum];
        uint16_t powerNumY = powerNumsY[testNum];
        double powerNumX2 = powerNumsX2[testNum];
        double powerNumY2 = powerNumsY2[testNum];
        size_t numThreads = 1;
        
        //x^powerNumX2 * cos(y^powerNumY)
        //y^powerNumY2 * cos(x^powerNumX)
        //Roots are at y^powerNumY = k1*pi + pi/2, x^powerNumX = k2*pi + pi/2
        //There are (upperBoundY^powerNumY)(upperBoundX^powerNumX)/pi^2 roots
        std::string functionString1 = "x^" + std::to_string(powerNumX2) + "*cos(y^" + std::to_string(powerNumY) +  ")"; //xcos(y)
        std::string functionString2 = "y^" + std::to_string(powerNumY2) + "*cos(x^" + std::to_string(powerNumX) +  ")"; //ycos(x)

        //Declare the intervals
        Interval startInterval;
        startInterval.lowerBounds.push_back(0.5);
        startInterval.lowerBounds.push_back(0.5);
        startInterval.upperBounds.push_back(upperBoundX);
        startInterval.upperBounds.push_back(upperBoundY);

        //Set up the solver
        std::vector<std::unique_ptr<Function>> functions;
        functions.emplace_back(std::make_unique<Function>(functionString1, variablesNames, subfunctionNames));
        functions.emplace_back(std::make_unique<Function>(functionString2, variablesNames, subfunctionNames));
        ThreadedSolver<Dimension::Two> solver(functions, numThreads, startInterval);
            
        //Solve it
        size_t trials = 1;
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        for(size_t i = 0; i < trials; i++) {
            solver.solve();
        }
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
        double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
        std::cout << "Solve " << functionString1 <<", "<< functionString2 << " on "<< startInterval.toString() << " with " << numThreads << " threads takes " <<nanos/(trials * 1e6)<< "ms.\n";
        
        //Assert the solutions
        std::vector<FoundRoot> foundRoots =  solver.getRoots();
        size_t expectedRootsX = power(upperBoundX, powerNumX) / M_PI + 0.5;
        size_t expectedRootsY = power(upperBoundX, powerNumX) / M_PI + 0.5;
        size_t expectedRoots = expectedRootsX * expectedRootsY;
        for(FoundRoot& root : foundRoots) {
            double rootModPiX = fmod(power(root.root[0], powerNumX), M_PI);
            double rootModPiY = fmod(power(root.root[1], powerNumY), M_PI);
            //std::cout<<rootModPiX-M_PI/2<<"\t"<<rootModPiY-M_PI/2<<"\n";
            XCTAssert(withinEpslion(rootModPiX, M_PI/2, 1e-4) && withinEpslion(rootModPiY, M_PI/2, 1e-4));
        }
        XCTAssert(expectedRoots == foundRoots.size());
        
        Timer::getTimingResultsAndClear();
    }
    std::cout<<"\n";
}

@end
