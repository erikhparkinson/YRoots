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

#include <iostream>
#include <fstream>
#include <vector>
#include "InputFileParser.h"
#include "ThreadedSolver.h"
#include "thread"
#include "Timer.h"

@interface TestThreadedSolver : XCTestCase

@end

@implementation TestThreadedSolver

- (void)setUp {
    //Disable the timer as default, as running with multiple threads is an issue
    Timer::disable();
    Function::clearSavedFunctions();
}

- (void)tearDown {
}

std::vector<std::vector<Function::SharedFunctionPtr>> createAllFunctions(const std::vector<std::string> functionStrings, const std::vector<std::string>& variablesNames, size_t numThreads) {
    Function::clearSavedFunctions();
    std::vector<std::vector<Function::SharedFunctionPtr>> result;
    result.resize(numThreads);
    for(size_t i = 0; i < functionStrings.size(); i++) {
        Function::addFunction("f" + std::to_string(i), functionStrings[i], variablesNames);
    }
    Function::addThreadFunctions(numThreads);
    
    for(size_t threadNum = 0; threadNum < numThreads; threadNum++) {
        for(size_t funcNum = 0; funcNum < functionStrings.size(); funcNum++) {
            result[threadNum].push_back(Function::getThreadFunctionByName(threadNum, "f" + std::to_string(funcNum)));
        }
    }

    return result;
}

- (void) testMain {
    Timer::enable();

    //Get the file names
    std::string inputFileName = "input.txt";
    
    InputFileParser inputParser (inputFileName);
    inputParser.parse();
        
    std::vector<std::vector<Function::SharedFunctionPtr>>& functions = inputParser.getFunctions();
    Interval interval = inputParser.getInterval();
    size_t numThreads = inputParser.getNumThreads();
    SubdivisionParameters subdivisionParameters = inputParser.getSubdivisionParameters();
    
    //Solve
    switch(functions[0].size()) {
        case 0: {
            printAndThrowRuntimeError("No functions found!");
            break;
        }
        case 1:
        {
            ThreadedSolver<Dimension::One> threadedSolver(functions, numThreads, interval, subdivisionParameters);
            threadedSolver.solve();
            break;
        }
        case 2:
        {
            ThreadedSolver<Dimension::Two> threadedSolver(functions, numThreads, interval, subdivisionParameters);
            threadedSolver.solve();
            break;
        }
        case 3:
        {
            ThreadedSolver<Dimension::Three> threadedSolver(functions, numThreads, interval, subdivisionParameters);
            threadedSolver.solve();
            break;
        }
        default:
        {
            ThreadedSolver<Dimension::NDim> threadedSolver(functions, numThreads, interval, subdivisionParameters);
            threadedSolver.solve();
            break;
        }
    }
    Timer::getTimingResultsAndClear();
}

- (void)testBasic1D {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    std::string functionString = "-1+2*x1^20";
    std::vector<std::string> functionStrings;
    functionStrings.push_back(functionString);
    
    Interval startInterval;
    startInterval.lowerBounds.push_back(-1.0);
    startInterval.upperBounds.push_back(1.0);
    
    SubdivisionParameters subdivisionParameters;
    
    std::cout<<"\n";
    for(size_t numThreads = 1; numThreads <= 4; numThreads++) {
        std::vector<std::vector<Function::SharedFunctionPtr>> functions = createAllFunctions(functionStrings, variablesNames, numThreads);
        ThreadedSolver<Dimension::One> solver(functions, numThreads, startInterval, subdivisionParameters);
            
        size_t trials = 1;
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        for(size_t i = 0; i < trials; i++) {
            solver.solve();
        }
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

        double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
        std::cout << "Solve with " << numThreads << " threads takes " << formatTimePretty(nanos/trials)<< ".\n";

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
    
    SubdivisionParameters subdivisionParameters;

    std::cout<<"\n";
    for(size_t testNum = 0; testNum < upperBounds.size(); testNum++) {
        //Get the variables
        double upperBound = upperBounds[testNum];
        double powerNum = powerNums[testNum];
        for(size_t numThreads = 1; numThreads <= 4; numThreads++) {
            std::vector<std::string> variablesNames;
            variablesNames.push_back("x");
            //sin(x^powerNum)
            //Roots at x^powerNum = k*pi
            //There are upperBound^powerNum/pi roots
            std::string functionString = "sin(x^" + std::to_string(powerNum) + ")";
            std::vector<std::string> functionStrings;
            functionStrings.push_back(functionString);
            
            //Declare the intervals
            Interval startInterval;
            startInterval.lowerBounds.push_back(0.5);
            startInterval.upperBounds.push_back(upperBound);

            //Set up the solver
            std::vector<std::vector<Function::SharedFunctionPtr>> functions = createAllFunctions(functionStrings, variablesNames, numThreads);
            ThreadedSolver<Dimension::One> solver(functions, numThreads, startInterval, subdivisionParameters);
                
            //Solve it
            size_t trials = 1;
            std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
            for(size_t i = 0; i < trials; i++) {
                solver.solve();
            }
            std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
            double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        
            std::cout << "Solve " << functionString << " on "<< startInterval.toString() << " with " << numThreads << " threads takes " << formatTimePretty(nanos/trials)<< ".\n";

            //Assert the solutions
            std::vector<FoundRoot> foundRoots =  solver.getRoots();
            size_t expectedRoots = power(upperBound, powerNum) / M_PI;
            for(FoundRoot& root : foundRoots) {
                double rootModPi = fmod(power(root.root[0], powerNum), M_PI);
                XCTAssert(withinEpslion(rootModPi, 0.0) || withinEpslion(rootModPi, M_PI));
            }
            XCTAssert(expectedRoots == foundRoots.size());
        }
    }
    std::cout<<"\n";
}

- (void)testHarder2D {
    //Parameters
    std::vector<double> upperBoundsX;
    std::vector<double> upperBoundsY;
    std::vector<double> powerNumsX;
    std::vector<double> powerNumsY;
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
    
    SubdivisionParameters subdivisionParameters;

    std::cout<<"\n";
    for(size_t testNum = 0; testNum < upperBoundsX.size(); testNum++) {
        //Get the variables
        double upperBoundX = upperBoundsX[testNum];
        double upperBoundY = upperBoundsY[testNum];
        double powerNumX = powerNumsX[testNum];
        double powerNumY = powerNumsY[testNum];
        double powerNumX2 = powerNumsX2[testNum];
        double powerNumY2 = powerNumsY2[testNum];
        for(size_t numThreads = 1; numThreads <= 4; numThreads++) {
            //x^powerNumX2 * cos(y^powerNumY)
            //y^powerNumY2 * cos(x^powerNumX)
            //Roots are at y^powerNumY = k1*pi + pi/2, x^powerNumX = k2*pi + pi/2
            //There are (upperBoundY^powerNumY)(upperBoundX^powerNumX)/pi^2 roots
            std::string functionString1 = "x^" + std::to_string(powerNumX2) + "*cos(y^" + std::to_string(powerNumY) +  ")"; //xcos(y)
            std::string functionString2 = "y^" + std::to_string(powerNumY2) + "*cos(x^" + std::to_string(powerNumX) +  ")"; //ycos(x)
            std::vector<std::string> functionStrings;
            functionStrings.push_back(functionString1);
            functionStrings.push_back(functionString2);
            
            //Declare the intervals
            Interval startInterval;
            startInterval.lowerBounds.push_back(0.5);
            startInterval.lowerBounds.push_back(0.5);
            startInterval.upperBounds.push_back(upperBoundX);
            startInterval.upperBounds.push_back(upperBoundY);

            //Set up the solver
            std::vector<std::vector<Function::SharedFunctionPtr>> functions = createAllFunctions(functionStrings, variablesNames, numThreads);
            ThreadedSolver<Dimension::Two> solver(functions, numThreads, startInterval, subdivisionParameters);
                
            //Solve it
            size_t trials = 1;
            std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
            for(size_t i = 0; i < trials; i++) {
                solver.solve();
            }
            std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
            double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        
            std::cout << "Solve " << functionString1 <<", "<< functionString2 << " on "<< startInterval.toString() << " with " << numThreads << " threads takes " << formatTimePretty(nanos/trials)<< ".\n";
            
            //Assert the solutions
            std::vector<FoundRoot> foundRoots =  solver.getRoots();
            size_t expectedRootsX = power(upperBoundX, powerNumX) / M_PI + 0.5;
            size_t expectedRootsY = power(upperBoundY, powerNumY) / M_PI + 0.5;
            size_t expectedRoots = expectedRootsX * expectedRootsY;
            for(FoundRoot& root : foundRoots) {
                double rootModPiX = fmod(power(root.root[0], powerNumX), M_PI);
                double rootModPiY = fmod(power(root.root[1], powerNumY), M_PI);
                XCTAssert(withinEpslion(rootModPiX, M_PI/2, 1e-4) && withinEpslion(rootModPiY, M_PI/2, 1e-4));
            }
            XCTAssert(expectedRoots == foundRoots.size());
        }
    }
    std::cout<<"\n";
}

- (void)testCustom {
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x");
    variablesNames.push_back("y");

    std::string functionString1 = "sin(30*x-y/30)+y";
    std::string functionString2 = "cos(x/30-30*y)-x";
    std::vector<std::string> functionStrings;
    functionStrings.push_back(functionString1);
    functionStrings.push_back(functionString2);
    
    //Declare the intervals
    Interval startInterval;
    startInterval.lowerBounds.push_back(-1);
    startInterval.lowerBounds.push_back(-1);
    startInterval.upperBounds.push_back(1);
    startInterval.upperBounds.push_back(1);
    
    SubdivisionParameters subdivisionParameters;

    size_t numThreads = 1;
    
    //Set up the solver
    std::vector<std::vector<Function::SharedFunctionPtr>> functions = createAllFunctions(functionStrings, variablesNames, numThreads);
    ThreadedSolver<Dimension::Two> solver(functions, numThreads, startInterval, subdivisionParameters);
    
    //Solve it
    size_t trials = 1;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        solver.solve();
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    std::cout << "Solve " << functionString1 <<", "<< functionString2 << " on "<< startInterval.toString() << " with " << numThreads << " threads takes " <<formatTimePretty(nanos/trials)<< ".\n";
    
    //Assert the solutions
    std::vector<FoundRoot> foundRoots =  solver.getRoots();
    XCTAssert(363 == foundRoots.size());
}

@end
