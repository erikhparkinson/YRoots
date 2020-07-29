//
//  ThreadedSolverND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ThreadedSolverND_ipp
#define ThreadedSolverND_ipp

template <Dimension D>
ThreadedSolver<D>::ThreadedSolver(const std::vector<std::vector<std::unique_ptr<FunctionInterface>>>& _allFunctions, size_t _numThreads, Interval& startInterval) :
m_allFunctions(_allFunctions),
m_numThreads(_numThreads),
m_intervalTracker(false),
m_rootTracker()
{
    //Limited the number of threads by the hardware amount
    m_numThreads = std::min(m_numThreads, (size_t)std::thread::hardware_concurrency());
    m_numRunningThreads.store(0);
    
    //TODO: Parse these parameters
    SubdivisionParameters subdivisionParameters;
    
    //Push back the first interval
    SolveParameters firstSolveParameters(startInterval, 0, std::vector<size_t>(_allFunctions[0].size(), subdivisionParameters.approximationDegree));
    m_firstSolveParameters = firstSolveParameters;
    //m_intervalsToRun.push(firstSolveParameters);
    
    m_killThreads.store(false);
    for(int threadNum = 0; threadNum < m_numThreads - 1; threadNum++) {
        //Create the solver
        m_subdivisionSolvers.emplace_back(std::make_unique<SubdivisionSolver<D>>(m_allFunctions[threadNum], m_intervalsToRun, subdivisionParameters, m_intervalTracker, m_rootTracker));

        //Create the threads
        //m_threadPool.push_back(std::thread(&ThreadedSolver<D>::runThread, this, threadNum));
    }
    
    //Run the main thread as well
    int threadNum = static_cast<int>(m_numThreads - 1);
    //Create the solver
    m_subdivisionSolvers.emplace_back(std::make_unique<SubdivisionSolver<D>>(m_allFunctions[threadNum], m_intervalsToRun, subdivisionParameters, m_intervalTracker, m_rootTracker));

    //Create the threads
    //runThread(threadNum);
    
    //Wait for the other threads to join
    for(std::thread &thread : m_threadPool) {
        thread.join();
    }
    m_threadPool.clear();
    
    //TODO: Have these print to a file
    //m_rootTracker.printResults();
    //m_intervalTracker.printResults();
}

template <Dimension D>
ThreadedSolver<D>::~ThreadedSolver() {
    //Kill the threads
    m_killThreads.store(true);
    //Wait for them to finish
    for(std::thread &thread : m_threadPool) {
        thread.join();
    }
    m_threadPool.clear();
}

template <Dimension D>
void ThreadedSolver<D>::solve() {
    m_intervalsToRun.push(m_firstSolveParameters);
    for(int threadNum = 0; threadNum < m_numThreads - 1; threadNum++) {
        //Create the threads
        m_threadPool.push_back(std::thread(&ThreadedSolver<D>::runThread, this, threadNum));
    }
    int threadNum = static_cast<int>(m_numThreads - 1);
    runThread(threadNum);
    for(std::thread &thread : m_threadPool) {
        thread.join();
    }
    m_threadPool.clear();
    
    m_rootTracker.printResults();
    m_intervalTracker.printResults();
}


template <Dimension D>
void ThreadedSolver<D>::runThread(size_t threadNum){
    //End the loop when kill threads is flipped
    while(!m_killThreads.load()) {
        //Try to get a parameter
        SolveParameters nextParameters;
        if(m_intervalsToRun.try_pop(nextParameters)) {
            m_numRunningThreads.fetch_add(1);
            m_subdivisionSolvers[threadNum]->solve(nextParameters);
            m_numRunningThreads.fetch_sub(1);
        }
        else if(m_numRunningThreads.load() == 0) {
            break;
        }
    }
}


#endif /* ThreadedSolverND_ipp */
