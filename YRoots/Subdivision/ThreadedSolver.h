//
//  ThreadedSolver.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ThreadedSolver_h
#define ThreadedSolver_h

#include <thread>
#include <atomic>
#include "../Subdivision/SubdivisionSolver.h"
#include "../SolutionTracking/IntervalTracker.h"
#include "../SolutionTracking/RootTracker.h"
#include "../Utilities/MultiPool.h"
#include "../Utilities/ConcurrentStack.h"

template <Dimension D>
class ThreadedSolver {
public:
    ThreadedSolver(std::vector<std::vector<Function::SharedFunctionPtr>>& _functions, const GeneralParameters& _generalParameters, Interval& _startInterval, const SubdivisionParameters& _subdivisionParameters);
    
    void solve();
    
    std::vector<FoundRoot> getRoots() {
        return m_rootTracker.getRoots();
    }
    
    ~ThreadedSolver();
    
private:
    void runThread(size_t _threadNum);
    
private:
    std::vector<std::vector<Function::SharedFunctionPtr>>   m_allFunctions;
    size_t                                                  m_numThreads;
    std::atomic<bool>                                       m_killThreads;
    std::atomic<size_t>                                     m_numRunningThreads;
    std::vector<std::thread>                                m_threadPool;
    
    //Thread safe data:
    ConcurrentStack<SolveParameters>                        m_intervalsToRun;
    MultiPool<SolveParameters>                              m_solveParametersPool;
    
    //For Storing Data
    IntervalTracker                                         m_intervalTracker;
    RootTracker                                             m_rootTracker;
    
    
    std::vector<std::unique_ptr<SubdivisionSolver<D>>>      m_subdivisionSolvers;
    
    SolveParameters*                                        m_firstSolveParameters;
};

#include "ThreadedSolverND.ipp"

#endif /* ThreadedSolver_h */
