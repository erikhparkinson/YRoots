//
//  TestQuadraticCheck.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 11/30/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "TestUtils.hpp"
#include "IntervalChecking/IntervalChecker.hpp"

template <Dimension D>
class IntervalCheckerMock : public IntervalChecker<D> {
public:
    IntervalCheckerMock(size_t _rank, IntervalTracker& _intervalTracker, size_t _threadNum, ConcurrentStack<SolveParameters>& _intervalsToRun, ObjectPool<SolveParameters>& _solveParametersPool) :
    IntervalChecker<D>(_rank, _intervalTracker, _threadNum, _intervalsToRun, _solveParametersPool)
    {

    }
    
    void runQuadraticCheck(ChebyshevApproximation<D>& _approximation) {
        IntervalChecker<D>::runQuadraticCheck(_approximation);
    }
    
    std::vector<bool>& get_m_intervalMask() {
        return IntervalChecker<D>::m_intervalMask;
    }
};


@interface TestQuadraticCheck : XCTestCase

@end

@implementation TestQuadraticCheck

- (void)setUp {
    m_allocated = false;
    Function::clearSavedFunctions();
}

- (void)tearDown {
    if(m_allocated) {
        //Deallocate everything
        fftw_free(m_approximation);
        m_allocated = false;
    }
}

void allocateMemoryTestQuadraticCheck()
{
    if(m_allocated) {
        //Deallocate everything
        fftw_free(m_approximation);
        m_allocated = false;
    }
    
    m_allocated = true;
    
    //Allocate the memory for the test once m_rank and m_approximationDegree have been set
    m_sideLength = 2*(m_approximationDegree+1);
    m_arrayLength = power(m_sideLength, m_rank);
    
    //Alllocate memory
    m_approximation = fftw_alloc_real(m_arrayLength);
}

- (void)test1DBasic {
    //Set up the Mock Class
    size_t rank = 1;
    IntervalTracker intervalTracker(rank, 0, false, false, 4.0);
    size_t threadNum = 0;
    ConcurrentStack<SolveParameters> intervalsToRun(1);
    SolveParameters defualtParams;
    ObjectPool<SolveParameters> solveParametersPool(defualtParams, 32);
    IntervalCheckerMock<Dimension::One> intervalChecker(rank, intervalTracker, threadNum, intervalsToRun, solveParametersPool);

    //Set up the Chebyshev approximation
    m_rank = 1;
    m_approximationDegree = 10;
    allocateMemoryTestQuadraticCheck();
    ChebyshevApproximation<Dimension::One> chebApproximation(m_rank);
    chebApproximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);

    //Initialize the coefficients as all ones
    for(size_t i = 0; i <= m_approximationDegree; i++) {
        m_approximation[i] = 1.0;
    }
    m_approximation[0] = 100.0;

    intervalChecker.runQuadraticCheck(chebApproximation);
    std::vector<bool>& intervalMask = intervalChecker.get_m_intervalMask();
    XCTAssert(intervalMask[0]);
    XCTAssert(intervalMask[1]);
}

- (void)testQuadCheckTiming {
    //Set up the Mock Class
    size_t rank = 1;
    IntervalTracker intervalTracker(rank, 0, false, false, 4.0);
    size_t threadNum = 0;
    ConcurrentStack<SolveParameters> intervalsToRun(1);
    SolveParameters defualtParams;
    ObjectPool<SolveParameters> solveParametersPool(defualtParams, 32);
    IntervalCheckerMock<Dimension::One> intervalChecker(rank, intervalTracker, threadNum, intervalsToRun, solveParametersPool);

    //Set up the Chebyshev approximation
    m_rank = 1;
    m_approximationDegree = 10;
    allocateMemoryTestQuadraticCheck();
    ChebyshevApproximation<Dimension::One> chebApproximation(m_rank);
    chebApproximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);

    //Initialize the coefficients as all ones
    for(size_t i = 0; i <= m_approximationDegree; i++) {
        m_approximation[i] = 1;
    }
    m_approximation[0] = 100;

    size_t trials = 1000000;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        intervalChecker.runQuadraticCheck(chebApproximation);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    std::cout << "\nQuadratic Check takes " << formatTimePretty(nanos/trials) << ".\n\n";

    std::vector<bool>& intervalMask = intervalChecker.get_m_intervalMask();
    XCTAssert(intervalMask[0]);
    XCTAssert(intervalMask[1]);
}


@end
