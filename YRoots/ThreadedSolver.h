//
//  ThreadedSolver.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#include <tbb/concurrent_queue.h>
#include <thread>
#include <atomic>
#include "SubdivisionSolver.h"
#include "IntervalTracker.h"
#include "RootTracker.h"

#ifndef ThreadedSolver_h
#define ThreadedSolver_h

template <Dimension D>
class ThreadedSolver {
public:
    ThreadedSolver(const std::vector<std::vector<std::unique_ptr<FunctionInterface>>>& _allFunctions, size_t _numThreads, Interval& startInterval);
    
    void solve();
    
    ~ThreadedSolver();
    
private:
    void runThread(size_t threadNum);
    
private:
    const std::vector<std::vector<std::unique_ptr<FunctionInterface>>>&   m_allFunctions;
    size_t                                                  m_numThreads;
    std::atomic<bool>                                       m_killThreads;
    std::atomic<size_t>                                     m_numRunningThreads;
    std::vector<std::thread>                                m_threadPool;
    
    //TODO: Make a buffer pool of SolveParameters and make intervals to run
    //a concurrent_queue of pointers to it. That will avoid a lot of copying.
    //The buffer pool can allocate 100 SolveParameters at first and make a
    //concurrent_queue of the ones available, [0-99]. When one is needed it gets
    //popped and then pushed back in when done. If nothing is there when popped
    //then allocate another 100 and push those numbers into the queue.
    tbb::strict_ppl::concurrent_queue<SolveParameters>      m_intervalsToRun;
    IntervalTracker                                         m_intervalTracker;
    RootTracker                                             m_rootTracker;
    
    tbb::concurrent_vector<std::unique_ptr<SubdivisionSolver<D>>>      m_subdivisionSolvers;
    
    SolveParameters                                                    m_firstSolveParameters;
};

#include "ThreadedSolverND.ipp"

#endif /* ThreadedSolver_h */
