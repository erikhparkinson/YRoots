//
//  IntervalData.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/6/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalData_h
#define IntervalData_h

#include "utilities.h"
#include <fstream>
#include <numeric>

struct IntervalResult {
    Interval    m_interval;
    SolveMethod m_howFound;
    
    IntervalResult() {}
    IntervalResult(Interval _interval, SolveMethod _howFound):
    m_interval(_interval),
    m_howFound(_howFound)
    {}
};

class IntervalTracker {
public:
    IntervalTracker(size_t _rank, size_t _numThreads, bool _enabled, double totalArea):
    m_rank(_rank),
    m_numThreads(_numThreads),
    m_enabled(_enabled),
    m_totalArea(totalArea),
    m_unitIntervalArea(power(2, _rank)),
    m_updatingProgressBar(0),
    m_lastAreaSolved(0),
    m_printedTooDeepWarning(false)
    {
        m_intervalResults.resize(m_numThreads);
        m_outputFile = "intervals.txt";
        m_areaSolved.resize(m_numThreads, 0);
    }
    
    void storeResult(size_t _threadNum, Interval& _interval, SolveMethod _howFound, double newSize = 0.0) {
        //TODO: Have this be two seperate levels of enabled
        if (m_enabled) {
            //New Size is on the [-1,1] scale. So take (1-newSize/2^n) * _interval.getArea()
            double areaSolved = (1 - newSize / m_unitIntervalArea) * _interval.getArea();
            assert(areaSolved > 0);
            m_areaSolved[_threadNum] += areaSolved;
            //Also, maybe have an option that counts intervals but doesn't store them.
            //I could even template the IntervalTracker of the trackingLevel to avoid checking m_enabled each time
            m_intervalResults[_threadNum].emplace_back(_interval, _howFound);
#ifndef TESTING
            updateProgressBar();
#endif
            
            if(unlikely(_howFound == SolveMethod::TooDeep && !m_printedTooDeepWarning.exchange(true))) {
                printWarning("Max Depth recusrion depth reached on solve, program may never finish. Try running with higher tolerances!");
            }
        }
    }
    
    void printResults() {
        for(size_t i = 0; i < m_numThreads; i++) {
            std::cout<< "Thread " << i << " solved " << m_intervalResults[i].size() << " intervals\n";
        }
        return;
    }
    
    void printResults2() {
        std::vector<size_t> reasons;
        for(size_t i = 0; i < 10; i++) {
            reasons.push_back(0);
        }
        
        for(size_t i = 0; i < m_numThreads; i++) {
            for(size_t j = 0; j < m_intervalResults[i].size(); j++) {
                reasons[m_intervalResults[i][j].m_howFound]++;
            }
        }
        
        for(size_t i = 0; i < 10; i++) {
            std::cout<<"Reason "<<i<<" count "<<reasons[i]<<"\n";
        }
    }
    
    void logResults() {
        std::ofstream file;
        int precision = std::numeric_limits<double>::digits10 + 1;
        file.open (m_outputFile);
        for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
            file<<"Thread " << threadNum << " solved " << m_intervalResults[threadNum].size() << " intervals.\n";
            for(size_t intervalNum = 0; intervalNum < m_intervalResults[threadNum].size(); intervalNum++) {
                //Log the Lower Bound
                file<<"[";
                for(size_t i = 0; i < m_intervalResults[threadNum][intervalNum].m_interval.lowerBounds.size(); i++) {
                    file<<std::setprecision(precision)<<m_intervalResults[threadNum][intervalNum].m_interval.lowerBounds[i];
                    if(i + 1 < m_intervalResults[threadNum][intervalNum].m_interval.lowerBounds.size()){
                        file<<",";
                    }
                }
                file<<"],";
                //Log the Upper Bound
                file<<"[";
                for(size_t i = 0; i < m_intervalResults[threadNum][intervalNum].m_interval.upperBounds.size(); i++) {
                    file<<std::setprecision(precision)<<m_intervalResults[threadNum][intervalNum].m_interval.upperBounds[i];
                    if(i + 1 < m_intervalResults[threadNum][intervalNum].m_interval.upperBounds.size()){
                        file<<",";
                    }
                }
                file<<"]\t";
                //Log the reason
                file<<m_intervalResults[threadNum][intervalNum].m_howFound<<"\n";
            }
        }
        file.close();
    }

    
private:
    void updateProgressBar() {
        if(m_updatingProgressBar.fetch_add(1) == 0) {
            int currentArea = std::round(100 * std::accumulate(m_areaSolved.begin(), m_areaSolved.end(), 0.0) / m_totalArea);
            assert(currentArea >= 0);
            assert(currentArea <= 100);
            if(currentArea > m_lastAreaSolved) {
                //Print the progress
                int barWidth = 70;
                std::cout << "\r[";
                int pos = barWidth * currentArea/100;
                for (int i = 0; i < barWidth; ++i) {
                    if (i < pos) std::cout << "=";
                    else if (i == pos) std::cout << ">";
                    else std::cout << " ";
                }
                std::cout << "] " << currentArea << " %";
                if(unlikely(currentArea == 100)) {
                    std::cout<<"\n"; //Print a newline at the end so the program ends on a newline.
                }
                std::cout.flush();
                
                m_lastAreaSolved = currentArea;
            }
        }
        m_updatingProgressBar.fetch_add(-1);
    }
    
    size_t                                      m_rank;
    size_t                                      m_numThreads;
    bool                                        m_enabled;
    std::vector<std::vector<IntervalResult>>    m_intervalResults;
    std::string                                 m_outputFile;
    double                                      m_totalArea;
    double                                      m_unitIntervalArea;
    std::vector<double>                         m_areaSolved;
    int                                         m_lastAreaSolved;
    std::atomic<int>                            m_updatingProgressBar;
    std::atomic<bool>                           m_printedTooDeepWarning;
};

#endif /* IntervalData_h */
