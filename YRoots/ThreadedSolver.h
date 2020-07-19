//
//  ThreadedSolver.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#include "tbb/concurrent_queue.h"
#include <thread>
#include <atomic>
#include "SubdivisionSolver.h"

#ifndef ThreadedSolver_h
#define ThreadedSolver_h

template <Dimension D>
class ThreadedSolver {
public:
    ThreadedSolver(const std::vector<std::vector<std::unique_ptr<FunctionInterface>>>& _allFunctions, size_t _numThreads, Interval& startInterval);
    
    ~ThreadedSolver();
    
private:
    void runThread(size_t threadNum);
    
private:
    const std::vector<std::vector<std::unique_ptr<FunctionInterface>>>&   m_allFunctions;
    size_t                                                  m_numThreads;
    std::atomic<bool>                                       m_killThreads;
    std::atomic<size_t>                                     m_numRunningThreads;
    std::vector<std::thread>                                m_threadPool;
    
    tbb::strict_ppl::concurrent_queue<SolveParameters>      m_intervalsToRun;
    
    std::vector<std::unique_ptr<SubdivisionSolver<D>>>      m_subdivisionSolvers;
};

#include "ThreadedSolverND.ipp"

#endif /* ThreadedSolver_h */
