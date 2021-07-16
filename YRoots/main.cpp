//
//  main.cpp
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#include "Solve.hpp"
#include "Utilities/Timer.hpp"

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

    //Call the solver
    mainSolver(inputFileName);
    
#ifdef USE_TIMING
    Timer::getTimingResultsAndClear();
#endif
    
    return 0;
}
