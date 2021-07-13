//
//  TestFunction.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 10/14/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "TestUtils.hpp"
#include "Functions/Function.hpp"

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
    variableNames.push_back("x");
    
    auto evalChecker1D = [&](const std::string& functionString, double evalPoint, double correct) {
        Function tempFunction("", functionString, variableNames);
        std::vector<double> inputPoints;
        inputPoints.push_back(evalPoint);
        double result = tempFunction.evaluate(inputPoints);
        if(!withinEpslion(result, correct)) {
            std::cout<<"Fail: " << functionString << "\t" <<result << "\t" << correct << "\n";
        }
        XCTAssert(withinEpslion(result, correct));
    };
    
    evalChecker1D("5+2*x", 5, 15);
    evalChecker1D("5+8+13", 5, 26);
    evalChecker1D("sin(x)*cos(x)", 5, sin(5)*cos(5));
    evalChecker1D("5*4*3*cosh(x)", 5, 60*cosh(5));
    evalChecker1D("5/8/log(x)", 5, .625/log(5));
    evalChecker1D("5*x^2", 5, 125);
    evalChecker1D("5*x**2", 5, 125);
    evalChecker1D("2*T3(x)", 2, 52);
    evalChecker1D("2*T3(0)", 2, 0); //TODO: Make sure this is interpreted as a constant.
    evalChecker1D("-3*T12(x)", cos(12), -3*cos(144));
    evalChecker1D("2*x^4*x^5", 2, 1024);
    evalChecker1D("2*(x-1)", 5, 8);
    evalChecker1D("10^(-9)", 5, 1e-9);
    evalChecker1D("1-x^2", 2, -3);
    evalChecker1D("1-2^x", 2, -3);
    evalChecker1D("-cos(x)", 0, -1);
    //Test Scientific Notation
    evalChecker1D("1e-4*cos(x)", 0, 1e-4);
    evalChecker1D("1e4*cos(x)", 0, 1e4);
    evalChecker1D("1.2e-3-3", 5, (1.2e-3)-3);
    evalChecker1D("1.2e3-3", 5, (1.2e3)-3);
    evalChecker1D("1.e-3-3", 5, (1.e-3)-3);
    evalChecker1D("1.e3-3", 5, (1.e3)-3);
    evalChecker1D("-1.e-3-3", 5, (-1.e-3)-3);
    //Test splitting in parenthesis
    evalChecker1D("cos(x^2)", 3, cos(9));
    evalChecker1D("cos(x**2)", 3, cos(9));
    evalChecker1D("cos(x+2)", 3, cos(5));
    evalChecker1D("cos(x-2)", 3, cos(1));
    evalChecker1D("cos(x/2)", 3, cos(1.5));
    evalChecker1D("cos(x*2)", 3, cos(6));
    //Test logs
    evalChecker1D("log(x)", 3, log(3));
    evalChecker1D("log2(x)", 3, log2(3));
    evalChecker1D("log10(x)", 3, log10(3));
    //Test sum notation
    evalChecker1D("sum(sin(x+i),i,-5,5)", 5, sin(0)+sin(1)+sin(2)+sin(3)+sin(4)+sin(5)+sin(6)+sin(7)+sin(8)+sin(9)+sin(10));
    evalChecker1D("prod(sin(x+i),i,-4,2)", 5, sin(1)*sin(2)*sin(3)*sin(4)*sin(5)*sin(6)*sin(7));

    //TODO: Add tests that it fails to parse.
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
            std::vector<std::vector<double> > grid;
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
    std::vector<std::vector<double> > grid;
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
