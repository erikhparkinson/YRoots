//
//  TestThreadedSolver.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 7/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include "TestUtils.hpp"
#include "Subdivision/ThreadedSolver.hpp"
#include "IO/InputFileParser.hpp"
#include "Utilities/Timer.hpp"
#include "Utilities/utilities.hpp"
#include "Solve.hpp"

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

std::vector<std::vector<Function::SharedFunctionPtr> > createAllFunctions(const std::vector<std::string> functionStrings, const std::vector<std::string>& variablesNames, size_t numThreads) {
    Function::clearSavedFunctions();
    std::vector<std::vector<Function::SharedFunctionPtr> > result;
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
    //return;
    Timer::enable();

    //Get the file names
    std::string inputFileName = "input.txt";
    mainSolver(inputFileName);
    
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
    GeneralParameters generalParameters;

    std::cout<<"\n";
    for(size_t numThreads = 1; numThreads <= 4; numThreads++) {
        generalParameters.numThreads = numThreads;
        std::vector<std::vector<Function::SharedFunctionPtr> > functions = createAllFunctions(functionStrings, variablesNames, numThreads);
        ThreadedSolver<1> solver(functions, generalParameters, startInterval, subdivisionParameters);
            
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
    GeneralParameters generalParameters;

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
            generalParameters.numThreads = numThreads;
            std::vector<std::vector<Function::SharedFunctionPtr> > functions = createAllFunctions(functionStrings, variablesNames, numThreads);
            ThreadedSolver<1> solver(functions, generalParameters, startInterval, subdivisionParameters);
                
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
                XCTAssert(withinEpslion(rootModPi, 0.0, 1.e-5) || withinEpslion(rootModPi, M_PI, 1.e-5));
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
    GeneralParameters generalParameters;

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
            if(numThreads == 1) {
                Timer::enable();
            }
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
            generalParameters.numThreads = numThreads;
            std::vector<std::vector<Function::SharedFunctionPtr> > functions = createAllFunctions(functionStrings, variablesNames, numThreads);
            ThreadedSolver<2> solver(functions, generalParameters, startInterval, subdivisionParameters);
                
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
                XCTAssert(withinEpslion(rootModPiX, M_PI/2, 1e-8) && withinEpslion(rootModPiY, M_PI/2, 1e-8));
            }
            XCTAssert(expectedRoots == foundRoots.size());
            
            if(numThreads == 1) {
                Timer::getTimingResultsAndClear();
            }
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
    GeneralParameters generalParameters;

    generalParameters.numThreads = 1;
    
    //Set up the solver
    std::vector<std::vector<Function::SharedFunctionPtr> > functions = createAllFunctions(functionStrings, variablesNames, generalParameters.numThreads);
    ThreadedSolver<2> solver(functions, generalParameters, startInterval, subdivisionParameters);
    
    //Solve it
    size_t trials = 1;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        solver.solve();
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    std::cout << "Solve " << functionString1 <<", "<< functionString2 << " on "<< startInterval.toString() << " with " << generalParameters.numThreads << " threads takes " <<formatTimePretty(nanos/trials)<< ".\n";
    
    //Assert the solutions
    std::vector<FoundRoot> foundRoots =  solver.getRoots();
    XCTAssert(363 == foundRoots.size());
}

std::vector<std::vector<double> > parseRootsFile(const std::string& _yrootsFile) {
    std::vector<std::vector<double> > roots;
    
    std::ifstream inputFile;
    std::string line;
    inputFile.open(_yrootsFile);
    while(std::getline(inputFile, line)) {
        std::vector<double> root;
        std::vector<std::string> rootString = split(line, ",");
        for(auto&& s : rootString) {
            root.push_back(std::stod(s));
        }
        roots.push_back(root);
    }
    inputFile.close();

    return roots;
}

double distanceBetweenPoints(const std::vector<double>& p1, const std::vector<double>& p2) {
    assert(p1.size() == p2.size());
    
    double val = 0;
    for(size_t i = 0; i < p1.size(); i++) {
        val += power(p1[i] - p2[i], 2);
    }
    return sqrt(val);
}

bool compareTestFiles(const std::string& _yrootsFile, const std::string& _chebRootsFile, std::vector<Function::SharedFunctionPtr>& functions) {
    //const double maxDistanceAllowed = 1e-5;
    const double maxResidualAllowed = 1e-5;

    std::vector<std::vector<double> > myRoots = parseRootsFile(_yrootsFile);
    std::vector<std::vector<double> > chebRoots = parseRootsFile(_chebRootsFile);
    
    if(myRoots.size() != chebRoots.size()) {
        std::cout<<"Unequal number of roots found!\n";
        return false;
    }
    const size_t numRoots = myRoots.size();
    
    std::set<double> spotsFound;
    std::vector<double> distances(numRoots);
    std::vector<std::vector<double> > yrootsResiduals(numRoots);
    std::vector<std::vector<double> > chebResiduals(numRoots);
    for(size_t rootSpot1 = 0; rootSpot1 < numRoots; rootSpot1++) {
        const std::vector<double>& myRoot = myRoots[rootSpot1];
        //Find the closest root
        double bestDistance = std::numeric_limits<double>::max();
        size_t bestIndex = -1;
        for(size_t rootSpot2 = 0; rootSpot2 < numRoots; rootSpot2++) {
            double d = distanceBetweenPoints(myRoot, chebRoots[rootSpot2]);
            if(d < bestDistance) {
                bestDistance = d;
                bestIndex = rootSpot2;
            }
        }
        distances[rootSpot1] = bestDistance;
        //Track it
        if(spotsFound.find(bestIndex) != spotsFound.end()) {
            std::cout<<"Couldn't match up roots!!\n";
            return false;
        }
        else {
            spotsFound.insert(bestIndex);
        }
        //Get the residuals
        const std::vector<double>& chebRoot = chebRoots[bestIndex];
        for(size_t funcNum = 0; funcNum < functions.size(); funcNum++) {
            yrootsResiduals[rootSpot1].push_back(std::abs(functions[funcNum]->evaluate<double>(myRoot)));
            chebResiduals[rootSpot1].push_back(std::abs(functions[funcNum]->evaluate<double>(chebRoot)));
        }
    }
    
    //Check the distances and residuals
    double maxResidual = 0;
    for(size_t rootSpot = 0; rootSpot < numRoots; rootSpot++) {
        /*if(distances[rootSpot] > maxDistanceAllowed) {
            return false;
        }*/
        for(size_t i = 0; i < functions.size(); i++) {
            maxResidual = std::max(maxResidual, yrootsResiduals[rootSpot][i]);
            if(yrootsResiduals[rootSpot][i] > maxResidualAllowed) {
                return false;
            }
        }

        /*std::cout << "D1: " << distances[rootSpot] << "\n";
        std::cout << "R1: ";
        for(size_t i = 0; i < functions.size(); i++) {
            std::cout << yrootsResiduals[rootSpot][i] << "\t";
        }
        std::cout<<"\n";
        std::cout << "R2: ";
        for(size_t i = 0; i < functions.size(); i++) {
            std::cout << chebResiduals[rootSpot][i] << "\t";
        }
        std::cout<<"\n";*/
    }
    
    std::cout<<"\tPass with max residual: "<<maxResidual<<"\n";
    
    return true;
}

- (void)testChebSuite {
    Timer::enable();
    //Read through all the ChebSuite files, run the test, and then check the results against what ChebFun gets.
    const std::string testFolder = "TestFiles/ChebSuiteTests/";
    const std::string testResultFolder = "TestFiles/ChebSuiteResults/";
    std::vector<std::string> testNames;
    testNames.push_back("test_1.1");
    testNames.push_back("test_1.2");
    testNames.push_back("test_1.3");
    testNames.push_back("test_1.4");
    testNames.push_back("test_1.5");
    testNames.push_back("test_2.1");
    testNames.push_back("test_2.2");
    testNames.push_back("test_2.3");
    testNames.push_back("test_2.4");
    testNames.push_back("test_2.5");
    testNames.push_back("test_3.1");
    testNames.push_back("test_3.2");
    testNames.push_back("test_4.1");
    testNames.push_back("test_4.2");
    testNames.push_back("test_5.1");
    //testNames.push_back("test_6.1"); //This one has a double root. //TODO: Find a way to make it pass if it finds a double root at 0 or one root at 0.
    testNames.push_back("test_6.2");
    testNames.push_back("test_6.3");
    testNames.push_back("test_7.1");
    testNames.push_back("test_7.2");
    testNames.push_back("test_7.3");
    testNames.push_back("test_7.4");
    testNames.push_back("test_8.1");
    testNames.push_back("test_8.2");
    testNames.push_back("test_9.1");
    testNames.push_back("test_9.2");
    testNames.push_back("test_10.1");

    for(size_t testNum = 0; testNum < testNames.size(); testNum++) {
        const std::string testName = testNames[testNum];
        std::cout << "Running Test " + testName << "\n";

        //Get file names
        std::vector<std::string> splitName = split(testName, "_");
        const std::string testInputFile = testFolder + testName + ".txt";
        const std::string testResultFile = testResultFolder + splitName[0] + "_roots_" + splitName[1] + ".csv";
        const std::string testOutputFile = "roots.csv";

        //Solve it
        Function::clearSavedFunctions();
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        mainSolver(testInputFile);
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
        double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        std::cout << "Solves in " <<formatTimePretty(nanos)<< ".\n";
        
        //Grab functions
        Function::clearSavedFunctions();
        InputFileParser inputParser(testInputFile);
        inputParser.parse();
        std::vector<Function::SharedFunctionPtr>& functions = inputParser.getFunctions()[0];

        //Check it
        XCTAssert(compareTestFiles(testOutputFile, testResultFile, functions));
    }
    Timer::getTimingResultsAndClear();
}

- (void)testDemoNotebook {
    Timer::enable();
    //Read through all the DemoNotebook examples, run the test, and then check the results against what our python code gets.
    const std::string testFolder = "TestFiles/DemoNotebookTests/";
    const std::string testResultFolder = "TestFiles/DemoNotebookResults/";
    std::vector<std::string> testNames;
    testNames.push_back("test_1");
    testNames.push_back("test_2");
    testNames.push_back("test_3");
    testNames.push_back("test_Devil");
    testNames.push_back("test_Rosenbrock");
    testNames.push_back("test_Trefethen");
    testNames.push_back("test_4");
    //testNames.push_back("test_5"); //Add this in once I have ND quad check or Lipshitz stuff.

    for(size_t testNum = 0; testNum < testNames.size(); testNum++) {
        const std::string testName = testNames[testNum];
        std::cout << "Running Test " + testName << "\n";

        //Get file names
        std::vector<std::string> splitName = split(testName, "_");
        const std::string testInputFile = testFolder + testName + ".txt";
        const std::string testResultFile = testResultFolder + splitName[0] + "_roots_" + splitName[1] + ".csv";
        const std::string testOutputFile = "roots.csv";

        //Solve it
        Function::clearSavedFunctions();
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        mainSolver(testInputFile);
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
        double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        std::cout << "Solves in " <<formatTimePretty(nanos)<< ".\n";
        
        //Grab functions
        Function::clearSavedFunctions();
        InputFileParser inputParser(testInputFile);
        inputParser.parse();
        std::vector<Function::SharedFunctionPtr>& functions = inputParser.getFunctions()[0];

        //Check it
        XCTAssert(compareTestFiles(testOutputFile, testResultFile, functions));
    }
    Timer::getTimingResultsAndClear();
}

- (void)testRandomStuffTiming {
    return;
    //Vars to run the test on
    std::vector<double> approximation(10000,0.0);
    size_t rank = 2;
    size_t sideLength = 20;
    size_t partialSideLength = 10;
    double sumAbsVal = 0.0;
    
    //Solve it
    size_t trials = 1;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t trailNum = 0; trailNum < trials; trailNum++) {
        //The code to test

        //Set up the needed variables
        std::vector<size_t> inputSpot(rank,0);
        std::vector<size_t> multipliers(rank, 1);
        for(size_t i = 1; i < rank; i++) {
            multipliers[i] = sideLength*multipliers[i-1];
        }

        //Iterate through all the combinations
        size_t spotToInc = 0;
        sumAbsVal += std::abs(approximation[0]);
        while (spotToInc < rank) {
            bool firstPass = true;
            while(++inputSpot[spotToInc] < partialSideLength) {
                size_t spot = 0;
                for (size_t i = 0; i < rank; i++) {
                    spot += inputSpot[i]*multipliers[i];//TODO: get rid of this for loop, inc the spot as we go. Replace this slow piece of code all over
                }
                std::cout<<"Spot " << spot << "\n";
                sumAbsVal += std::abs(approximation[spot]);
                
                if(firstPass && spotToInc != 0) {
                    spotToInc = 0;
                }
                firstPass = false;
            }
            inputSpot[spotToInc] = 0;
            spotToInc++;
        }
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    std::cout << "This code takes " << formatTimePretty(nanos/trials)<< ".\n";
}


@end
