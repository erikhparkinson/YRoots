//
//  Solve.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/24/21.
//  Copyright © 2021 Erik Hales Parkinson. All rights reserved.
//

#ifndef Solve_h
#define Solve_h

#include <iostream>
#include <fstream>
#include <vector>
#include <thread>

#include "IO/InputFileParser.h"
#include "Subdivision/ThreadedSolver.h"
#include "Functions/Function.h"
#include "Utilities/Timer.h"

template<Dimension D>
void runSolve(std::vector<std::vector<Function::SharedFunctionPtr> >& _functions, const GeneralParameters& _generalParameters, Interval& _interval, const SubdivisionParameters& _subdivisionParameters) {
    Timer& m_timer = Timer::getInstance();
    size_t mainConstructorIndex = -1;
    size_t mainSolveIndex = -1;
    m_timer.registerTimer(mainConstructorIndex, "Main Constructor");
    m_timer.registerTimer(mainSolveIndex, "Main Solve");

    m_timer.startTimer(mainConstructorIndex);
    ThreadedSolver<D> threadedSolver(_functions, _generalParameters, _interval, _subdivisionParameters);
    m_timer.stopTimer(mainConstructorIndex);
    m_timer.startTimer(mainSolveIndex);
    threadedSolver.solve();
    m_timer.stopTimer(mainSolveIndex);
}

void mainSolver(const std::string& inputFileName) {    
    InputFileParser inputParser(inputFileName);
    inputParser.parse();
        
    std::vector<std::vector<Function::SharedFunctionPtr> >& functions = inputParser.getFunctions();
    Interval& interval = inputParser.getInterval();
    const GeneralParameters& generalParameters = inputParser.getGeneralParameters();
    const SubdivisionParameters& subdivisionParameters = inputParser.getSubdivisionParameters();

    #ifdef USE_TIMING
    if(!inputParser.getGeneralParameters().useTimer) {
        Timer::disable();
    }
    #endif
    
    //Solve
    switch(functions[0].size()) {
        case 0: {
            printAndThrowRuntimeError("No functions found!");
            break;
        }
        case 1:
        {
            runSolve<Dimension::One>(functions, generalParameters, interval, subdivisionParameters);
            break;
        }
        case 2:
        {
            runSolve<Dimension::Two>(functions, generalParameters, interval, subdivisionParameters);
            break;
        }
        case 3:
        {
            runSolve<Dimension::Three>(functions, generalParameters, interval, subdivisionParameters);
            break;
        }
        default:
        {
            runSolve<Dimension::NDim>(functions, generalParameters, interval, subdivisionParameters);
            break;
        }
    }

}

#endif /* Solve_h */
