//
//  testPolynomial.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 9/2/21.
//  Copyright © 2021 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "TestUtils.hpp"
#include "Functions/Function.hpp"
#include "Utilities/ErrorTracker.hpp"
#include <stdlib.h>
#include <time.h>

@interface testPolynomial : XCTestCase

@end

@implementation testPolynomial

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testPoly1 {
    Polynomial myPoly;
    myPoly.setRank(1);
    std::vector<Monomial> monomials;
    Monomial monomial;
    monomial.spot.resize(1);
    monomial.spot[0] = 0; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 1; monomial.coeff = 2;
    monomials.push_back(monomial);
    monomial.spot[0] = 3; monomial.coeff = 3;
    monomials.push_back(monomial);
    
    myPoly.addMonomials(monomials);
    myPoly.prepEvaluation();

    std::vector<double> evalPoints;
    evalPoints.push_back(5);
    
    double eval = myPoly.evaluate<double>(evalPoints);
    XCTAssert(withinEpslion(eval, 386));
}

- (void)testPoly2 {
    Polynomial myPoly;
    myPoly.setRank(2);
    std::vector<Monomial> monomials;
    Monomial monomial;
    monomial.spot.resize(2);
    monomial.spot[0] = 0; monomial.spot[1] = 1; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 1; monomial.spot[1] = 1; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 3; monomial.spot[1] = 1; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 2; monomial.spot[1] = 3; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 3; monomial.spot[1] = 3; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 4; monomial.spot[1] = 3; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 0; monomial.spot[1] = 4; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 1; monomial.spot[1] = 4; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 4; monomial.spot[1] = 4; monomial.coeff = 1;
    monomials.push_back(monomial);
    
    myPoly.addMonomials(monomials);
    myPoly.prepEvaluation();

    std::vector<double> evalPoints;
    evalPoints.push_back(2);
    evalPoints.push_back(2);

    double eval = myPoly.evaluate<double>(evalPoints);
    XCTAssert(withinEpslion(eval, 550));
}

- (void)testPoly3 {
    Polynomial myPoly;
    myPoly.setRank(2);
    std::vector<Monomial> monomials;
    Monomial monomial;
    monomial.spot.resize(2);
    monomial.spot[0] = 1; monomial.spot[1] = 0; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 1; monomial.spot[1] = 1; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 1; monomial.spot[1] = 3; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 3; monomial.spot[1] = 2; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 3; monomial.spot[1] = 3; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 3; monomial.spot[1] = 4; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 4; monomial.spot[1] = 0; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 4; monomial.spot[1] = 1; monomial.coeff = 1;
    monomials.push_back(monomial);
    monomial.spot[0] = 4; monomial.spot[1] = 4; monomial.coeff = 1;
    monomials.push_back(monomial);
    myPoly.addMonomials(monomials);
    myPoly.prepEvaluation();
    
    std::vector<double> evalPoints;
    evalPoints.push_back(2);
    evalPoints.push_back(2);

    double eval = myPoly.evaluate<double>(evalPoints);
    XCTAssert(withinEpslion(eval, 550));
}

double randomUniform1() { //Random uniform [0,1].
    return static_cast<double>(rand()) / RAND_MAX;
}

double randomUniform2() { //Random uniform [-1,1].
    return 2 * randomUniform1() - 1;
}

- (void)testPoly2DHard {
    srand(1923717); //Seed the randomness
    
    const size_t rank = 2;
    
    //The Testing Parameters
    std::vector<size_t> polyDegreesToTest;
    polyDegreesToTest.push_back(2);
    polyDegreesToTest.push_back(3);
    polyDegreesToTest.push_back(5);
    polyDegreesToTest.push_back(10);
    polyDegreesToTest.push_back(25);

    std::vector<size_t> gridSizesToTest;
    gridSizesToTest.push_back(2);
    gridSizesToTest.push_back(5);
    gridSizesToTest.push_back(10);

    std::vector<double> densitiesToTest;
    densitiesToTest.push_back(0.01);
    densitiesToTest.push_back(0.1);
    densitiesToTest.push_back(0.5);
    densitiesToTest.push_back(0.9);
    densitiesToTest.push_back(0.99);
    densitiesToTest.push_back(1.0);
    
    const size_t loopsPerTest = 100;
    const double epsilon = 1e-5;

    for(size_t& polyDegree : polyDegreesToTest) {
        for(size_t& gridSize : gridSizesToTest) {
            for(double& density : densitiesToTest) {
                std::cout<<"Testing Poly Evals Degree: " << polyDegree << ", Grid Size: " << gridSize << ", Density: " << density << "\n";
                for(size_t testNum = 0; testNum < loopsPerTest; testNum++) {
                    //std::cout<<"Testing Loop: " << testNum << "\n";
                    //Create the objects
                    Polynomial myPoly;
                    myPoly.setRank(rank);
                    Monomial myMonomial;
                    myMonomial.spot.resize(rank);
                    
                    //Create a random polynomial
                    for(size_t i = 0; i < polyDegree; i++) {
                        myMonomial.spot[0] = i;
                        for(size_t j = 0; j < polyDegree; j++) {
                            myMonomial.spot[1] = j;
                            myMonomial.coeff = randomUniform2(); //Random uniform [-1,1].
                            if(randomUniform1() <= density) {
                                myPoly.addMonomial(myMonomial);
                            }
                        }
                    }
                    myPoly.prepEvaluation();
                    
                    //Create a random grid
                    std::vector<std::vector<double>> myGrid;
                    for(size_t i = 0; i < rank; i++) {
                        myGrid.resize(i+1);
                        for(size_t j = 0; j < gridSize; j++) {
                            myGrid[i].push_back(randomUniform2()); //Random uniform [-1,1].
                        }
                    }
                    
                    //Evaluate the grid
                    std::vector<double> gridResults(power(gridSize, rank), 0);
                    myPoly.evaluateGrid(myGrid, gridResults);
                    
                    //Evaluate over every point seperately
                    std::vector<double> evalPoint(rank, 0);
                    size_t gridPoint = 0;
                    if(myPoly.getNumUsedDimensions() == 1) {
                        size_t spotToEval = myPoly.getHasDimension()[0] ? 0 : 1;
                        for(size_t j = 0; j < gridSize; j++) {
                            evalPoint[spotToEval] = myGrid[spotToEval][j];
                            double result = myPoly.evaluate<double>(evalPoint); //Quick Eval
                            double result2 = myPoly.evaluateSlow(evalPoint); //Slow Eval to compare
                            double result3 = gridResults[gridPoint++]; //Grid Eval
                            //std::cout << result << "\t" << result2 << "\t"<< result3 << "\n";
                            XCTAssert(withinEpslion(result, result2, epsilon));
                            XCTAssert(withinEpslion(result, result3, epsilon));
                        }
                    }
                    else {
                        for(size_t i = 0; i < rank; i++) {
                            evalPoint[1] = myGrid[1][i];
                            for(size_t j = 0; j < gridSize; j++) {
                                evalPoint[0] = myGrid[0][j];
                                double result = myPoly.evaluate<double>(evalPoint); //Quick Eval
                                double result2 = myPoly.evaluateSlow(evalPoint); //Slow Eval to compare
                                double result3 = gridResults[gridPoint++]; //Grid Eval
                                //std::cout << result << "\t" << result2 << "\t"<< result3 << "\n";
                                XCTAssert(withinEpslion(result, result2, epsilon));
                                XCTAssert(withinEpslion(result, result3, epsilon));
                            }
                        }
                    }
                }
            }
        }
    }
}

//TODO: Add a 3D Test

- (void)testPoly2DTiming {
    const size_t polyDegree = 100;
    const bool triangular = false;
    const size_t gridSize = 10;
    const double density = 1.0;
    const size_t rank = 2;
    
    //Create the objects
    Polynomial myPoly;
    myPoly.setRank(rank);
    Monomial myMonomial;
    myMonomial.spot.resize(rank);
    
    //Create a random polynomial
    for(size_t i = 0; i < polyDegree; i++) {
        myMonomial.spot[0] = i;
        for(size_t j = 0; j < polyDegree; j++) {
            myMonomial.spot[1] = j;
            myMonomial.coeff = randomUniform2(); //Random uniform [-1,1].
            if(randomUniform1() <= density) {
                if(i + j <= polyDegree || !triangular) {
                    myPoly.addMonomial(myMonomial);
                }
            }
        }
    }
    size_t prepTrials = 1;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < prepTrials; i++) {
        myPoly.prepEvaluation();
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    std::cout << "\nEvaluation Prep takes for Dimension 2 Degree " << polyDegree << " Density " << density << " takes " << formatTimePretty(nanos/prepTrials)<< ".\n";
    
    //Create a random grid
    std::vector<std::vector<double>> myGrid;
    for(size_t i = 0; i < rank; i++) {
        myGrid.resize(i+1);
        for(size_t j = 0; j < gridSize; j++) {
            myGrid[i].push_back(randomUniform2()); //Random uniform [-1,1].
        }
    }
    
    //Evaluate the grid
    std::vector<double> gridResults(power(gridSize, rank), 0);
    
    size_t evalTrials = 1000;
    start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < evalTrials; i++) {
        myPoly.evaluateGrid(myGrid, gridResults);
    }
    end = std::chrono::high_resolution_clock::now();
    nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    std::cout << "Evaluate Grid of size " << gridSize << " takes for Dimension 2 Degree " << polyDegree << " Density " << density << " takes " << formatTimePretty(nanos/evalTrials)<< ".\n\n";
}
    
//TODO: Add tests that test random polys.
//Have the 2D test, for a bunch of different powers and density probabilities, create random polys with randomly filled spots
//Evaluate the monomials in a simple for loop on a grid, and then call evaluate grid and evaluate on each point.

@end
