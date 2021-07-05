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
    Function::clearSavedFunctions();
}

- (void)tearDown {

}

- (void)testFunction1DBasic {
    std::vector<double> inputPoints;
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    
    std::string functionString = "5+2*x0";
    Function tempFunction("", functionString, variableNames);
    inputPoints.push_back(5);
    double result = tempFunction.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 15));
    
    functionString = "5+8+13";
    Function tempFunction2("", functionString, variableNames);
    result = tempFunction2.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 26));

    functionString = "sin(x0)*cos(x0)";
    Function tempFunction3("", functionString, variableNames);
    result = tempFunction3.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, sin(5)*cos(5)));

    functionString = "5*4*3*cosh(x0)";
    Function tempFunction4("", functionString, variableNames);
    result = tempFunction4.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 60*cosh(5)));
    
    functionString = "5/8/ln(x0)";
    Function tempFunction5("", functionString, variableNames);
    result = tempFunction5.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, .625*log(5)));
    
    functionString = "log(x0,125)";
    Function tempFunction6("", functionString, variableNames);
    result = tempFunction6.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 3.0));

    functionString = "log(5,x0)";
    Function tempFunction7("", functionString, variableNames);
    result = tempFunction7.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 1.0));

    functionString = "5*x0^2";
    Function tempFunction8("", functionString, variableNames);
    result = tempFunction8.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 125.));

    functionString = "2*T3(x0)";
    Function tempFunction9("", functionString, variableNames);
    inputPoints[0] = 2;
    result = tempFunction9.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 52));

    functionString = "-3*T12(x0)";
    Function tempFunction10("", functionString, variableNames);
    inputPoints[0] = cos(12);
    result = tempFunction10.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, -3*cos(144)));

    functionString = "2*x0^4*x0^5";
    Function tempFunction11("", functionString, variableNames);
    inputPoints[0] = 2;
    result = tempFunction11.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 1024));

    functionString = "2*(x0-1)";
    Function tempFunction12("", functionString, variableNames);
    inputPoints[0] = 5;
    result = tempFunction12.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 8));

    functionString = "10^(-9)";
    Function tempFunction13("", functionString, variableNames);
    inputPoints[0] = 5;
    result = tempFunction13.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 1e-9));

    functionString = "1-x0^2";
    Function tempFunction14("", functionString, variableNames);
    inputPoints[0] = 2;
    result = tempFunction14.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, -3));

    functionString = "1-2x0^2";
    Function tempFunction15("", functionString, variableNames);
    inputPoints[0] = 2;
    result = tempFunction15.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, -7));

    functionString = "1-2^x0";
    Function tempFunction16("", functionString, variableNames);
    inputPoints[0] = 2;
    result = tempFunction16.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, -3));
}

- (void)testFunction2DBasic {
    std::vector<double> inputPoints;
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");

    std::string functionString = "5+2*x0+x1";
    Function tempFunction("", functionString, variableNames);
    inputPoints.push_back(3);
    inputPoints.push_back(8);
    double result = tempFunction.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 19));

    functionString = "5+2*x0^3*x1^2";
    Function tempFunction2("", functionString, variableNames);
    result = tempFunction2.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 3461));

    functionString = "5+2*x0^2*x1^2*x0";
    Function tempFunction3("", functionString, variableNames);
    result = tempFunction3.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 3461));

    functionString = "5+2*T2(x0)*T2(x0)";
    Function tempFunction4("", functionString, variableNames);
    result = tempFunction4.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 583));

    functionString = "5+2*T2(x0)*(T2(x1)*T2(x0))";
    Function tempFunction5("", functionString, variableNames);
    result = tempFunction5.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 73411));

    functionString = "5+2*T2(x0)*(T2(x0)*T2(x1))";
    Function tempFunction6("", functionString, variableNames);
    result = tempFunction6.evaluate(inputPoints);
    XCTAssert(withinEpslion(result, 73411));
}

- (void)testEvalGrid2D{
    std::vector<std::string> variableNames;
    variableNames.push_back("x");
    variableNames.push_back("y");
    
    std::vector<std::string> functionStrings;
    functionStrings.push_back("7");
    functionStrings.push_back("x");
    functionStrings.push_back("y");
    functionStrings.push_back("x-y");
    functionStrings.push_back("x-y+.5");
    functionStrings.push_back("x+y");
    functionStrings.push_back("sin(30*x-y/30)+y");
    functionStrings.push_back("cos(x/30-30*y)-x");
    functionStrings.push_back("x^1.3*cos(y^0.8)");
    functionStrings.push_back("cos(x*y)/sin(x-14.8*y+e)-log(14)+(1+x+2*y)^tan(2)");
    functionStrings.push_back("1+x^40+y^40");

    for(size_t functionNumber = 0; functionNumber < functionStrings.size(); functionNumber++) {
        const std::string functionString = functionStrings[functionNumber];
        std::cout<<"Testing Grid Eval on " << functionString << "\n";

        Function tempFunction("", functionString, variableNames);
        
        size_t dimension = variableNames.size();
        for(size_t numPoints = 1; numPoints < 100; numPoints++) {
            std::vector<std::vector<double>> grid;
            grid.resize(dimension);
            for(size_t i = 0; i < dimension; i++) {
                for(size_t j = 0; j < numPoints; j++) {
                    grid[i].push_back(j);
                }
            }
            std::vector<double> results;
            results.resize(power(numPoints, dimension));

            tempFunction.evaluateGrid(grid, results);
            
            std::vector<double> evalPoints;
            evalPoints.resize(2);
            for(size_t j = 0; j < numPoints; j++) {
                for(size_t i = 0; i < numPoints; i++) {
                    evalPoints[0] = grid[0][i];
                    evalPoints[1] = grid[1][j];
                    double eval = tempFunction.evaluate(evalPoints);
                    if(!withinEpslion(eval, results[j*numPoints + i])) {
                        std::cout << eval << "\t" << results[j*numPoints + i] << "\t" << std::abs(eval-results[j*numPoints + i]) << "\n";
                    }
                    XCTAssert(withinEpslion(eval, results[j*numPoints + i]));
                }
            }
        }
    }
}

- (void)testFunctionTiming {
    for(size_t numFives = 1; numFives < 10; numFives++) {
        std::vector<double> inputPoints;
        std::vector<std::string> variableNames;
        variableNames.push_back("x0");
        std::string functionString = "5";
            for(size_t i = 1; i < numFives; i++) {
                functionString += "+5";
            }
        
        Function tempFunction("", functionString, variableNames);
        inputPoints.push_back(2);
        double result = 0;
        
        size_t trials = 1000000;
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        for(size_t i = 0; i < trials; i++) {
            result = tempFunction.evaluate(inputPoints);
        }
        std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

        double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        
        std::cout << "Evaluating function: " << numFives << " takes " << formatTimePretty(nanos/(trials)) << ".\n";
        //std::cout<<result<<"\n";
    }
}

- (void)testFunctionEvaluateGridTiming {
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");
    std::string functionString = "sin(cos(x0)+x1)-x0^2*x1^cos(tanh(x0/x1^3))+sin(cos(x0)+x1)+x0^2";
    
    size_t dimension = variableNames.size();
    size_t numPoints = 100;
    double numToEval = 1.7;
    
    Function tempFunction("", functionString, variableNames);
    std::vector<std::vector<double>> grid;
    grid.resize(dimension);
    for(size_t i = 0; i < dimension; i++) {
        for(size_t j = 0; j < numPoints; j++) {
            grid[i].push_back(numToEval);
        }
    }
    std::vector<double> results;
    results.resize(power(numPoints, dimension));
    
    size_t trials = 1000;
    tempFunction.evaluateGrid(grid, results);
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        tempFunction.evaluateGrid(grid, results);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "\nEvaluating function: " << functionString << " on grid takes " << formatTimePretty(nanos/trials)<< ".\n\n";
}

- (void)testFunctionParsingTiming {
    std::vector<std::string> variableNames;
    variableNames.push_back("x0");
    variableNames.push_back("x1");
    std::string functionString = "sin(cos(x0)+x1)-x0^2*x1^cos(tanh(x0/x1^3))-sin(cos(x0)+x1)+x0^2";
        
    size_t trials = 1000;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < trials; i++) {
        Function tempFunction("", functionString, variableNames);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "\nCreating function: " << functionString << " takes " << formatTimePretty(nanos/trials) << ".\n\n";
}



@end
