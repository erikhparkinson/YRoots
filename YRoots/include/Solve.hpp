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

#include "IO/InputFileParser.hpp"
#include "Subdivision/ThreadedSolver.hpp"
#include "Functions/Function.hpp"
#include "Utilities/Timer.hpp"

template<int Rank>
void runSolve(std::vector<std::vector<Function::SharedFunctionPtr> >& _functions, const GeneralParameters& _generalParameters, Interval& _interval, const SubdivisionParameters& _subdivisionParameters, const VariableSubsitutionInfo& _variableSubsitutionInfo) {
    Timer& m_timer = Timer::getInstance();
    static size_t mainConstructorIndex = -1;
    static size_t mainSolveIndex = -1;
    m_timer.registerTimer(mainConstructorIndex, "Main Constructor");
    m_timer.registerTimer(mainSolveIndex, "Main Solve");

    m_timer.startTimer(mainConstructorIndex);
    ThreadedSolver<Rank> threadedSolver(_functions, _generalParameters, _interval, _subdivisionParameters, _variableSubsitutionInfo);
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
    const VariableSubsitutionInfo& variableSubsitutionInfo = inputParser.getVariableSubstitutionInfo();

    #ifdef USE_TIMING
    if(!generalParameters.useTimer) {
        Timer::disable();
    }
    else {
        Timer::enable();
    }
    #endif
    
    //Solve
    switch(functions[0].size()) {
        case 0:
            printAndThrowRuntimeError("No functions found!");
            break;
        case 1:
            runSolve<1>(functions, generalParameters, interval, subdivisionParameters, variableSubsitutionInfo);
            break;
        case 2:
            runSolve<2>(functions, generalParameters, interval, subdivisionParameters, variableSubsitutionInfo);
            break;
        case 3:
            runSolve<3>(functions, generalParameters, interval, subdivisionParameters, variableSubsitutionInfo);
            break;
        default:
            runSolve<-1>(functions, generalParameters, interval, subdivisionParameters, variableSubsitutionInfo);
            break;
    }

}

#endif /* Solve_h */
