//
//  TestFunction.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 10/14/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "TestUtils.h"
#include "Function.h"

@interface TestFunction : XCTestCase

@end

@implementation TestFunction

- (void)setUp {

}

- (void)tearDown {

}

- (void)testFunction1DBasic {
    std::vector<double> inputPoints;
    std::vector<std::string> variableNames;
    std::vector<std::string> subfunctionNames;
    variableNames.push_back("x0");
    
    std::string functionString = "5+2*x0";
    Function tempFunction(functionString, variableNames, subfunctionNames);
    inputPoints.push_back(5);
    double result = tempFunction.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 15));
    
    functionString = "5+8+13";
    Function tempFunction2(functionString, variableNames, subfunctionNames);
    result = tempFunction2.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 26));

    functionString = "sin(x0)*cos(x0)";
    Function tempFunction3(functionString, variableNames, subfunctionNames);
    result = tempFunction3.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, sin(5)*cos(5)));

    functionString = "5*4*3*cosh(x0)";
    Function tempFunction4(functionString, variableNames, subfunctionNames);
    result = tempFunction4.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 60*cosh(5)));
    
    functionString = "5/8/ln(x0)";
    Function tempFunction5(functionString, variableNames, subfunctionNames);
    result = tempFunction5.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, .625*log(5)));
    
    functionString = "log(x0,125)";
    Function tempFunction6(functionString, variableNames, subfunctionNames);
    result = tempFunction6.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 3.0));

    functionString = "log(5,x0)";
    Function tempFunction7(functionString, variableNames, subfunctionNames);
    result = tempFunction7.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 1.0));

    functionString = "5*x0^2";
    Function tempFunction8(functionString, variableNames, subfunctionNames);
    result = tempFunction8.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 125.));

    functionString = "2*T3(x0)";
    Function tempFunction9(functionString, variableNames, subfunctionNames);
    inputPoints[0] = 2;
    result = tempFunction9.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 52));

    functionString = "-3*T12(x0)";
    Function tempFunction10(functionString, variableNames, subfunctionNames);
    inputPoints[0] = cos(12);
    result = tempFunction10.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, -3*cos(144)));

    functionString = "2*x0^4*x0^5";
    Function tempFunction11(functionString, variableNames, subfunctionNames);
    inputPoints[0] = 2;
    result = tempFunction11.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 1024));
}

- (void)testFunction2DBasic {
    std::vector<double> inputPoints;
    std::vector<std::string> variableNames;
    std::vector<std::string> subfunctionNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");

    std::string functionString = "5+2*x0+x1";
    Function tempFunction(functionString, variableNames, subfunctionNames);
    inputPoints.push_back(3);
    inputPoints.push_back(8);
    double result = tempFunction.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 19));

    functionString = "5+2*x0^3*x1^2";
    Function tempFunction2(functionString, variableNames, subfunctionNames);
    result = tempFunction2.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 3461));

    functionString = "5+2*x0^2*x1^2*x0";
    Function tempFunction3(functionString, variableNames, subfunctionNames);
    result = tempFunction3.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 3461));

    functionString = "5+2*T2(x0)*T2(x0)";
    Function tempFunction4(functionString, variableNames, subfunctionNames);
    result = tempFunction4.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 583));

    functionString = "5+2*T2(x0)*(T2(x1)*T2(x0))";
    Function tempFunction5(functionString, variableNames, subfunctionNames);
    result = tempFunction5.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 73411));

    functionString = "5+2*T2(x0)*(T2(x0)*T2(x1))";
    Function tempFunction6(functionString, variableNames, subfunctionNames);
    result = tempFunction6.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 73411));
}

- (void)testFunctionTiming {
    for(size_t numFives = 1; numFives < 10; numFives++) {
    std::vector<double> inputPoints;
    std::vector<std::string> variableNames;
    std::vector<std::string> subfunctionNames;
    variableNames.push_back("x0");
    std::string functionString = "5";
        for(size_t i = 1; i < numFives; i++) {
            functionString += "+5";
        }
    
    Function tempFunction(functionString, variableNames, subfunctionNames);
    inputPoints.push_back(2);
    double result = 0;
    
    size_t trials = 1000000;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        result = tempFunction.evaluate(inputPoints);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "Evaluating function: " << numFives << " takes " <<nanos/(trials)<< "ns.\n";
    //std::cout<<result<<"\n";
    }
}



@end
