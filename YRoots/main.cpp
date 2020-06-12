//
//  main.cpp
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "InputFileParser.h"
#include "SubdivisionSolver.h"

int main(int argc, const char * argv[]) {
    if(argc != 3) {
        std::cout<<"Arguments input and output must be given!\n";
        return 0;
    }
    //Get the file names
    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];
    
    InputFileParser inputParser (inputFileName);
    
    std::vector<std::unique_ptr<FunctionInterface>> functions = inputParser.parseFunctions();
    
    //TODO: Have an input for the search interval
    Interval interval;
    for(size_t i = 0; i < functions.size(); i++) {
        interval.lowerBounds.push_back(-1.0);
        interval.upperBounds.push_back(1.0);
    }
    
    //TODO: Check that the number of dimensions equals the number of variables
    switch(functions.size()) {
        case 0:
            throw std::runtime_error("No functions found!");
            break;
        case 1:
        {
            SubdivisionSolver<Dimension::One> subdivisionSolver(functions);
            subdivisionSolver.solve(interval, 0);
            break;
        }
        case 2:
        {
            SubdivisionSolver<Dimension::Two> subdivisionSolver(functions);
            subdivisionSolver.solve(interval, 0);
            break;
        }
        default:
        {
            SubdivisionSolver<Dimension::NDim> subdivisionSolver(functions);
            subdivisionSolver.solve(interval, 0);
            break;
        }
    }
    
    return 0;
}
