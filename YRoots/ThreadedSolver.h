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
#include "SubdivisionSolver.h"
#include "IntervalTracker.h"
#include "RootTracker.h"
#include "MultiPool.h"
#include "ConcurrentStack.h"

template <Dimension D>
class ThreadedSolver {
public:
    ThreadedSolver(const std::vector<std::vector<std::unique_ptr<FunctionInterface>>>& _allFunctions, size_t _numThreads, Interval& startInterval);
    
    void solve();
    
    ~ThreadedSolver();
    
private:
    void runThread(size_t _threadNum);
    
private:
    const std::vector<std::vector<std::unique_ptr<FunctionInterface>>>&   m_allFunctions;
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
