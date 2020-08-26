//
//  TestChebyshevApproximation.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 8/11/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "ChebyshevApproximation.hpp"

@interface TestChebyshevApproximation : XCTestCase

@end

@implementation TestChebyshevApproximation


- (void)setUp {
    m_allocated = false;
}

- (void)tearDown {
    if(m_allocated) {
        //Deallocate everything
        fftw_free(m_approximation);
        m_allocated = false;
    }
}

void allocateMemoryTestChebyshevApproximation()
{
    if(m_allocated) {
        //Deallocate everything
        fftw_free(m_approximation);
        m_allocated = false;
    }
    
    m_allocated = true;
    
    //Allocate the memory for the test once m_rank and m_approximationDegree have been set
    m_sideLength = 2*m_approximationDegree;
    m_arrayLength = power(m_sideLength, m_rank);
    
    //Alllocate memory
    m_approximation = fftw_alloc_real(m_arrayLength);
}

void fillApproximationTestChebyshevApproximation(double number) {
    for(size_t arraySpot = 0; arraySpot < m_arrayLength; arraySpot++) {
        bool fillSpot = true;
        size_t tempSpot = arraySpot;
        for(size_t dim = 0; dim < m_rank; dim++) {
            if(tempSpot % m_sideLength > m_approximationDegree) {
                fillSpot = false;
                break;
            }
            tempSpot /= m_sideLength;
        }
        if(fillSpot) {
            m_approximation[arraySpot] = number;
        }
    }
}

- (void) test1DTrimCoeff {
    m_rank = 1;
    m_approximationDegree = 10;
    allocateMemoryTestChebyshevApproximation();
    ChebyshevApproximation<Dimension::One> approximation(m_rank);
    
    //Initialize the coefficients as all ones
    for(size_t i = 0; i < m_arrayLength; i++) {
        m_approximation[i] = 1;
    }
    
    //We can approximate degree 10 as it's already that degree
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 10));
    //We can also do degree 100 and it's fine
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 100));

    //Test not being able to trim, and then being able to
    for(size_t targetDegree = m_approximationDegree - 1; targetDegree > 0; targetDegree--) {
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, targetDegree));
        m_approximation[targetDegree+1] = 1e-12;
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, targetDegree));
    }
    
    //Test trimming multiple degrees
    fillApproximationTestChebyshevApproximation(1e-12);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 1));

    //Test failing trimming multiple degrees
    fillApproximationTestChebyshevApproximation(4e-11);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, 1));
    
    //Test the initial error being too big
    fillApproximationTestChebyshevApproximation(1e-12);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 1e-9);
    XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, 1));

    //Test with relative tol and large inf norm
    fillApproximationTestChebyshevApproximation(1e-5);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1e6, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 1));

    //Test with large abs_tol
    fillApproximationTestChebyshevApproximation(1e-5);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-4, 1e-10, 1));
}

- (void) test2DTrimCoeff {
    m_rank = 2;
    m_approximationDegree = 9;
    allocateMemoryTestChebyshevApproximation();
    ChebyshevApproximation<Dimension::Two> approximation(m_rank);
    
    //Initialize the coefficients as all ones
    for(size_t i = 0; i < m_arrayLength; i++) {
        m_approximation[i] = 1;
    }
    
    //We can approximate degree 20 as it's already that degree
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 20));
    //We can also do degree 100 and it's fine
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 100));

    //Test not being able to trim, and then being able to
    for(size_t targetDegree = m_approximationDegree - 1; targetDegree > 0; targetDegree--) {
        ///TODO: Write this generically
        //approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        //XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, targetDegree));
        //m_approximation[targetDegree+1] = 1e-12;
        //approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        //XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, targetDegree));
    }
    
    //Test trimming multiple degrees
    fillApproximationTestChebyshevApproximation(1e-12);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 1));

    //Test failing trimming multiple degrees
    fillApproximationTestChebyshevApproximation(1e-11);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, 1));
    
    //Test the initial error being too big
    fillApproximationTestChebyshevApproximation(1e-12);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 1e-9);
    XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, 1));

    //Test with relative tol and large inf norm
    fillApproximationTestChebyshevApproximation(1e-5);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1e7, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 1));

    //Test with large abs_tol
    fillApproximationTestChebyshevApproximation(1e-5);
    approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
    XCTAssert(approximation.trimCoefficients(1e-3, 1e-10, 1));
}

- (void) testNDTrimCoeff {
    for(m_rank = 1; m_rank < 5; m_rank++) {
        m_approximationDegree = 9;
        allocateMemoryTestChebyshevApproximation();
        ChebyshevApproximation<Dimension::NDim> approximation(m_rank);
        
        double testMultipler = power(10.0, m_rank-1);
        
        //Initialize the coefficients as all ones
        for(size_t i = 0; i < m_arrayLength; i++) {
            m_approximation[i] = 1;
        }
        
        //We can approximate degree 20 as it's already that degree
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 10*m_rank));
        //We can also do degree 100 and it's fine
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 100*m_rank));
        
        //Test trimming multiple degrees
        fillApproximationTestChebyshevApproximation(1e-12/testMultipler);
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 1));

        //Test failing trimming multiple degrees
        fillApproximationTestChebyshevApproximation(1e-10/testMultipler);
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, 1));
        
        //Test the initial error being too big
        fillApproximationTestChebyshevApproximation(1e-12/testMultipler);
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 1e-9);
        XCTAssert(!approximation.trimCoefficients(1e-10, 1e-10, 1));

        //Test with relative tol and large inf norm
        fillApproximationTestChebyshevApproximation(1e-5);
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1e6*testMultipler, false, 0);
        XCTAssert(approximation.trimCoefficients(1e-10, 1e-10, 1));

        //Test with large abs_tol
        fillApproximationTestChebyshevApproximation(1e-5);
        approximation.setApproximation(m_approximationDegree, m_sideLength, m_approximation, 1, false, 0);
        XCTAssert(approximation.trimCoefficients(1e-4*testMultipler, 1e-10, 1));
    }
}


@end
