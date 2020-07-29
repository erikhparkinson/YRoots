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
#include "ThreadedSolver.h"
#include "thread"

int main(int argc, const char * argv[]) {
    if(argc != 3) {
        std::cout<<"Arguments input and output must be given!\n";
        return 0;
    }
    //Get the file names
    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];
    
    InputFileParser inputParser (inputFileName);
    
    //TODO: Have an input for number of threads
    size_t numThreads = 1;
    
    std::vector<std::vector<std::unique_ptr<FunctionInterface>>> functions = inputParser.parseFunctions(numThreads);
    
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
            ThreadedSolver<Dimension::One> threadedSolver(functions, numThreads, interval);
            threadedSolver.solve();
            break;
        }
        case 2:
        {
            ThreadedSolver<Dimension::Two> threadedSolver(functions, numThreads, interval);
            threadedSolver.solve();
            break;
        }
        default:
        {
            ThreadedSolver<Dimension::NDim> threadedSolver(functions, numThreads, interval);
            threadedSolver.solve();
            break;
        }
    }
    
    return 0;
}
