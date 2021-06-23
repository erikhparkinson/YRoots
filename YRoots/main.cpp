//
//  main.cpp
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef USE_TIMING
#define USE_TIMING
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include "InputFileParser.h"
#include "ThreadedSolver.h"
#include "thread"
#include "Timer.h"

template<Dimension D>
void runSolve(std::vector<std::vector<Function::SharedFunctionPtr>>& _functions, size_t _numThreads, Interval& _interval, const SubdivisionParameters& _subdivisionParameters) {
    Timer& m_timer = Timer::getInstance();
    size_t mainConstructorIndex = -1;
    size_t mainSolveIndex = -1;
    m_timer.registerTimer(mainConstructorIndex, "Main Constructor");
    m_timer.registerTimer(mainSolveIndex, "Main Solve");

    m_timer.startTimer(mainConstructorIndex);
    ThreadedSolver<D> threadedSolver(_functions, _numThreads, _interval, _subdivisionParameters);
    m_timer.stopTimer(mainConstructorIndex);
    m_timer.startTimer(mainSolveIndex);
    threadedSolver.solve();
    m_timer.stopTimer(mainSolveIndex);
}

int main(int argc, const char * argv[]) {
    #ifdef USE_TIMING
        Timer::enable();
    #endif

    //Get the file names
    std::string inputFileName;
    if(argc == 1) {
        inputFileName = "input.txt";
    }
    else if (argc == 2) {
        inputFileName = argv[1];
    }
    else {
        printAndThrowRuntimeError("Only one input allowed!");
    }
    
    InputFileParser inputParser (inputFileName);
    inputParser.parse();
        
    std::vector<std::vector<Function::SharedFunctionPtr>>& functions = inputParser.getFunctions();
    Interval interval = inputParser.getInterval();
    size_t numThreads = inputParser.getNumThreads();
    SubdivisionParameters subdivisionParameters = inputParser.getSubdivisionParameters();

    //Solve
    switch(functions[0].size()) {
        case 0: {
            printAndThrowRuntimeError("No functions found!");
            break;
        }
        case 1:
        {
            runSolve<Dimension::One>(functions, numThreads, interval, subdivisionParameters);
            break;
        }
        case 2:
        {
            runSolve<Dimension::Two>(functions, numThreads, interval, subdivisionParameters);
            break;
        }
        case 3:
        {
            runSolve<Dimension::Three>(functions, numThreads, interval, subdivisionParameters);
            break;
        }
        default:
        {
            runSolve<Dimension::NDim>(functions, numThreads, interval, subdivisionParameters);
            break;
        }
    }
    
    
#ifdef USE_TIMING
    Timer::getTimingResultsAndClear();
#endif
    
    return 0;
}
