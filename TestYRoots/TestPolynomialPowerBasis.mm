//
//  TestYRoots.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 5/30/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef TestPolynomialPowerBasis_mm
#define TestPolynomialPowerBasis_mm

#import <XCTest/XCTest.h>
#include "PowerBasisPolynomial.h"
#include <iostream>

@interface TestPolynomialPowerBasis : XCTestCase

@end

@implementation TestPolynomialPowerBasis

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

 //TODO: Add tests with more complex monomials
 
template<class T1, class T2>
bool withinEpslion(T1 a, T2 b, double epsilon = 1.e-10) {
    return std::abs(a-b) < epsilon;
}

- (void)testBasic1D {
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    std::string polynomialString = "5+2*x0";
    PowerBasisPolynomial testPoly(polynomialString, variableNames);
    std::vector<double> points;
    
    points.push_back(0);
    XCTAssert(withinEpslion(testPoly.evaluate(points), 5));
    points[0] = 1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 7));
    points[0] = 10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 25));
}

- (void)testBasicHigherDegree1D {
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    std::string polynomialString = "5+2*x0-4*x0^3-x0^7+x0^8";
    PowerBasisPolynomial testPoly(polynomialString, variableNames);
    std::vector<double> points;
    
    points.push_back(0);
    XCTAssert(withinEpslion(testPoly.evaluate(points), 5));
    points[0] = 1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 3));
    points[0] = 10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 89996025));
}

- (void)testBasic2D {
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");
    std::string polynomialString = "5+2*x0+3*x1";
    PowerBasisPolynomial testPoly(polynomialString, variableNames);
    std::vector<double> points;
    
    points.push_back(0);
    points.push_back(0);
    XCTAssert(withinEpslion(testPoly.evaluate(points), 5));
    points[0] = 1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 7));
    points[0] = 10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 25));
    points[1] = 1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 28));
    points[1] = 10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 55));
}

- (void)testBasicHigherDegree2D {
    // Test higher x0 power.
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");
    std::string polynomialString = "5+2*x0^2+3*x1";
    PowerBasisPolynomial testPoly1(polynomialString, variableNames);
    std::vector<double> points;
    
    points.push_back(0);
    points.push_back(0);
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 5));
    points[0] = 1;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 7));
    points[0] = 10;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 205));
    points[1] = 1;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 208));
    points[1] = 10;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 235));

    //Test higher x1 power
    polynomialString = "5+2*x0+3*x1^2";
    PowerBasisPolynomial testPoly2(polynomialString, variableNames);
    points[0] = 0;
    points[1] = 0;
    XCTAssert(withinEpslion(testPoly2.evaluate(points), 5));
    points[0] = 1;
    XCTAssert(withinEpslion(testPoly2.evaluate(points), 7));
    points[0] = 10;
    XCTAssert(withinEpslion(testPoly2.evaluate(points), 25));
    points[1] = 1;
    XCTAssert(withinEpslion(testPoly2.evaluate(points), 28));
    points[1] = 10;
    XCTAssert(withinEpslion(testPoly2.evaluate(points), 325));
}

- (void)testComplex2D {
    // Test higher x0 power.
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");
    std::string polynomialString = "5+2*x0^2+3*x1-2*x0*x1";
    PowerBasisPolynomial testPoly1(polynomialString, variableNames);
    std::vector<double> points;
    
    points.push_back(0);
    points.push_back(0);
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 5));
    points[0] = 1;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 7));
    points[0] = 10;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 205));
    points[1] = 1;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 188));
    points[1] = 10;
    XCTAssert(withinEpslion(testPoly1.evaluate(points), 35));
}

- (void)testBasicND {
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");
    variableNames.push_back("x2");
    variableNames.push_back("x3");
    variableNames.push_back("x4");
    variableNames.push_back("x5");
    variableNames.push_back("x6");
    std::string polynomialString = "-10+5*x6+2*x0-3*x1-4*x2+.5*x5+.1*x4-x3";
    PowerBasisPolynomial testPoly(polynomialString, variableNames);
    std::vector<double> points(7,0);
    
    XCTAssert(withinEpslion(testPoly.evaluate(points), -10));
    points[0] = 1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -8));
    points[0] = 10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 10));
    points[1] = 1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), 7));
    points[1] = 10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -20));
    points[2] = 1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -24));
    points[2] = 10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -60));
    points[3] = -1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -59));
    points[3] = -10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -50));
    points[4] = -1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -50.1));
    points[4] = -10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -51));
    points[5] = -1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -51.5));
    points[5] = -10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -56));
    points[6] = -1;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -61));
    points[6] = -10;
    XCTAssert(withinEpslion(testPoly.evaluate(points), -106));
}


@end

#endif
