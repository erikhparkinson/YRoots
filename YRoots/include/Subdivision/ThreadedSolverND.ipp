//
//  ThreadedSolverND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/18/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ThreadedSolverND_ipp
#define ThreadedSolverND_ipp

template <int Rank>
ThreadedSolver<Rank>::ThreadedSolver(std::vector<std::vector<Function::SharedFunctionPtr>>& _functions, const GeneralParameters& _generalParameters, Interval& _startInterval, const SubdivisionParameters& _subdivisionParameters, const VariableSubsitutionInfo& _variableSubsitutionInfo) :
m_allFunctions(_functions),
m_numThreads(_generalParameters.numThreads),
m_killThreads(false),
m_numRunningThreads(0),
m_intervalsToRun(m_numThreads),
m_intervalTracker(_functions[0].size(), _generalParameters, _startInterval.getArea()),
m_rootTracker(m_numThreads, m_allFunctions, _generalParameters, _variableSubsitutionInfo)
{
    if(m_numThreads == 0) {
        printAndThrowRuntimeError("Can't run with 0 threads!");
    }
    else if(_functions.size() != m_numThreads) {
        printAndThrowRuntimeError("Not enough functions passed to solver!");
    }
    
    //Create the Solve Parameters Pool
    size_t rank = _functions[0].size();
    SolveParameters defaultSolveParameters;
    defaultSolveParameters.interval.lowerBounds.resize(rank);
    defaultSolveParameters.interval.upperBounds.resize(rank);
    defaultSolveParameters.goodDegrees.resize(rank);
    for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
        m_solveParametersPool.emplace_back(defaultSolveParameters, 1024);
    }
    
    //Create the first solve parameters
    m_firstSolveParameters = m_solveParametersPool[0].pop();
    m_firstSolveParameters->interval = _startInterval;
    std::fill(m_firstSolveParameters->goodDegrees.begin(), m_firstSolveParameters->goodDegrees.end(), _subdivisionParameters.approximationDegree);
    m_firstSolveParameters->currentLevel = 0;

    //Create the solvers
    for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
        m_subdivisionSolvers.emplace_back(::make_unique<SubdivisionSolver<Rank>>(threadNum, m_allFunctions[threadNum], m_intervalsToRun, m_solveParametersPool[threadNum], _subdivisionParameters, m_intervalTracker, m_rootTracker));
    }
}

template <int Rank>
ThreadedSolver<Rank>::~ThreadedSolver() {
    //Kill the threads
    m_killThreads.store(true);
    //Wait for them to finish
    for(std::thread &thread : m_threadPool) {
        thread.join();
    }
    m_threadPool.clear();
}

template <int Rank>
void ThreadedSolver<Rank>::solve() {
    //TODO: Create the threads in the constructor and have the threads help create the m_subdivisionSolvers?
    //Also, why do I have the sleep for a millisecond here???
    
    m_intervalsToRun.push(0, m_firstSolveParameters);
    for(size_t threadNum = 0; threadNum + 1 < m_numThreads; threadNum++) {
        //Create the threads
        m_threadPool.push_back(std::thread(&ThreadedSolver<Rank>::runThread, this, threadNum));
        std::this_thread::sleep_for(std::chrono::microseconds(1500));
    }
    
    //Run the front thread as well
    int threadNum = static_cast<int>(m_numThreads - 1);
    runThread(threadNum);
    
    //Wait for the threads to finish
    for(std::thread &thread : m_threadPool) {
        thread.join();
    }
    m_threadPool.clear();
    
    m_rootTracker.logResults();
    m_intervalTracker.logResults();
}

template <int Rank>
void ThreadedSolver<Rank>::runThread(size_t _threadNum){
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
