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
ThreadedSolver<D>::ThreadedSolver(std::vector<std::unique_ptr<Function>>& _functions, size_t _numThreads, Interval& _startInterval) :
m_numThreads(_numThreads),
m_intervalsToRun(m_numThreads),
m_intervalTracker(m_numThreads, true, _startInterval.getArea()), //TODO: Whether to store should be a parameter
m_rootTracker(m_numThreads)
{
    //Limited the number of threads by the hardware amount
    m_numThreads = std::min(m_numThreads, (size_t)std::thread::hardware_concurrency());
    m_numRunningThreads.store(0);

    //Make m_allFunctions
    m_allFunctions.resize(m_numThreads);
    for(size_t allFuncNum = 0; allFuncNum + 1 < m_numThreads; allFuncNum++) {
        for(size_t funcNum = 0; funcNum < _functions.size(); funcNum++) {
            m_allFunctions[allFuncNum].emplace_back(std::make_unique<Function>(*_functions[funcNum].get()));
        }
    }
    //This final function we just move the pointer over to avoid an unneeded copy
    for(size_t funcNum = 0; funcNum < _functions.size(); funcNum++) {
        m_allFunctions[m_allFunctions.size() - 1].emplace_back(_functions[funcNum].release());
    }

    
    //TODO: Parse these parameters
    SubdivisionParameters subdivisionParameters;
        
    //Create the subdivision solvers
    m_killThreads.store(false);
    
    //Create the Solve Parameters Pool
    size_t rank = _functions.size();
    SolveParameters defaultSolveParameters;
    defaultSolveParameters.interval.lowerBounds.resize(rank);
    defaultSolveParameters.interval.upperBounds.resize(rank);
    defaultSolveParameters.goodDegrees.resize(rank);
    for(int threadNum = 0; threadNum < m_numThreads; threadNum++) {
        m_solveParametersPool.emplace_back(defaultSolveParameters, 1024);
    }
    
    //Create the first solve parameters
    m_firstSolveParameters = m_solveParametersPool[0].pop();
    m_firstSolveParameters->interval = _startInterval;
    std::fill(m_firstSolveParameters->goodDegrees.begin(), m_firstSolveParameters->goodDegrees.end(), subdivisionParameters.approximationDegree);
    m_firstSolveParameters->currentLevel = 0;

    //Create the solvers
    for(int threadNum = 0; threadNum < m_numThreads; threadNum++) {
        m_subdivisionSolvers.emplace_back(std::make_unique<SubdivisionSolver<D>>(threadNum, m_allFunctions[threadNum], m_intervalsToRun, m_solveParametersPool[threadNum], subdivisionParameters, m_intervalTracker, m_rootTracker));
    }
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
    
    m_intervalsToRun.push(0, m_firstSolveParameters);
    for(int threadNum = 0; threadNum < m_numThreads - 1; threadNum++) {
        //Create the threads
        m_threadPool.push_back(std::thread(&ThreadedSolver<D>::runThread, this, threadNum));
    }
    
    //Run the front thread as well
    int threadNum = static_cast<int>(m_numThreads - 1);
    runThread(threadNum);
    
    //Wait for the threads to finish
    for(std::thread &thread : m_threadPool) {
        thread.join();
    }
    m_threadPool.clear();
    
    //m_rootTracker.logResults();
    //m_intervalTracker.logResults();
}

template <Dimension D>
void ThreadedSolver<D>::runThread(size_t _threadNum){
    //End the loop when kill threads is flipped
    SolveParameters* parameters;
    while(!m_killThreads.load()) {
        //Try to get a parameter
        if((parameters = m_intervalsToRun.pop(_threadNum))) {
            m_numRunningThreads.fetch_add(1);
            m_subdivisionSolvers[_threadNum]->solve(parameters);
            m_numRunningThreads.fetch_sub(1);
            m_solveParametersPool[_threadNum].push(parameters);
        }
        else if(m_numRunningThreads.load() == 0) {
            break;
        }
    }
}


#endif /* ThreadedSolverND_ipp */
