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
#include <chrono>
#include "TestUtils.hpp"
#include "Approximation/IntervalApproximator.hpp"
#include "Functions/PowerBasisPolynomial.hpp"

@interface TestIntervalApproximater : XCTestCase

@end

@implementation TestIntervalApproximater

- (void)setUp {
    m_allocated = false;
    Function::clearSavedFunctions();
}

- (void)tearDown {
    if(m_allocated) {
        //Deallocate everything
        fftw_free(m_inputPartial);
        free(m_kinds);
        fftw_free(m_output);
        fftw_free(m_input);
        m_allocated = false;
    }
}

void allocateMemoryTestIntervalApproximater()
{
    m_allocated = true;
    
    //Allocate the memory for the test once m_rank and m_approximationDegree have been set
    m_sideLength = 2*m_approximationDegree;
    m_arrayLength = power(m_sideLength, m_rank);
    m_partialSideLength = m_approximationDegree + 1;
    m_partialArrayLength = power(m_partialSideLength, m_rank);
    
    //Alllocate memory
    m_input = fftw_alloc_real(m_arrayLength);
    m_output= fftw_alloc_real(m_arrayLength);
    m_kinds = (fftw_r2r_kind*) malloc(m_rank * sizeof (fftw_r2r_kind));
    m_inputPartial = fftw_alloc_real(m_partialArrayLength);

    //Define the kinds
    for(size_t i = 0; i < m_rank; i++) {
        m_kinds[i] = FFTW_R2HC;
    }
}

- (void)testBasic1D {
    m_rank = 1;
    m_approximationDegree = 3;
    allocateMemoryTestIntervalApproximater();
    
    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    std::string functionString = "5+x1^2";
    
    Function::SharedFunctionPtr function = std::make_shared<Function>("", functionString, variablesNames);
    IntervalApproximator<1> intervalApproximator(m_rank, m_approximationDegree, m_input, m_output, m_kinds, m_partialArrayLength);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0);
    intervalApproximator.approximate(function, currentInterval, true);

    XCTAssert(withinEpslion(intervalApproximator.getOutput()[0], 5.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[1], 0.0));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[2], 0.5));
    XCTAssert(withinEpslion(intervalApproximator.getOutput()[3], 0.0));
}

- (void)testBasic2D {
    m_rank = 2;
    m_approximationDegree = 3;
    allocateMemoryTestIntervalApproximater();

    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    variablesNames.push_back("x2");
    std::string functionString = "5+x1^2+x2";
    
    Function::SharedFunctionPtr function = std::make_shared<Function>("", functionString, variablesNames);
    IntervalApproximator<2> intervalApproximator(m_rank, m_approximationDegree, m_input, m_output, m_kinds, m_partialArrayLength);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0);
    intervalApproximator.approximate(function, currentInterval, true);

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
    m_rank = 3;
    m_approximationDegree = 2;
    allocateMemoryTestIntervalApproximater();

    std::vector<std::string> variablesNames;
    variablesNames.push_back("x1");
    variablesNames.push_back("x2");
    variablesNames.push_back("x3");
    std::string functionString = "5+x1^2+x2+5*x1*x2*x3";
    
    Function::SharedFunctionPtr function = std::make_shared<Function>("", functionString, variablesNames);
    IntervalApproximator<2> intervalApproximator(m_rank, m_approximationDegree, m_input, m_output, m_kinds, m_partialArrayLength);
    Interval currentInterval;
    currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0); currentInterval.lowerBounds.push_back(-1.0);
    currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0); currentInterval.upperBounds.push_back(1.0);
    intervalApproximator.approximate(function, currentInterval, true);
    
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

- (void)testTimingTemp {
    m_rank = 2;
    m_approximationDegree = 5;
    allocateMemoryTestIntervalApproximater();
    
    size_t degreePoly = 40;

    std::vector<std::string> variablesNames;
    std::string functionString = "1+";
    for(size_t i = 0; i < m_rank; i++) {
        variablesNames.push_back("x" + std::to_string(i));
        functionString += "x" + std::to_string(i) + "^" + std::to_string(degreePoly);
        if (i+1 != m_rank) {
            functionString += "*";
        }
    }
    
    Function::SharedFunctionPtr function = std::make_shared<Function>("", functionString, variablesNames);
    IntervalApproximator<2> intervalApproximator(m_rank, m_approximationDegree, m_input, m_output, m_kinds, m_partialArrayLength);
    Interval currentInterval;
    for(size_t i = 0; i < m_rank; i++) {
        currentInterval.lowerBounds.push_back(-1.0);
        currentInterval.upperBounds.push_back(1.0);
    }

    size_t trials = 1000;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        intervalApproximator.approximate(function, currentInterval, true);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "\nApproximating a Degree " << degreePoly << " Dimension " << m_rank << " Polynomial to degree " << m_approximationDegree << " takes " <<nanos/(trials*1000)<< "us.\n\n";
}


@end

#endif
