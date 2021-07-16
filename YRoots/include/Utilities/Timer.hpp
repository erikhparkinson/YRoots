//
//  TimingDetails.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 6/22/21.
//  Copyright © 2021 Erik Hales Parkinson. All rights reserved.
//

#ifndef Timer_h
#define Timer_h

#include <chrono>
#include <fstream>
#include "Utilities/utilities.hpp"

class TimingDetails {
    //TODO: Have this hold a vector of the individual run times, so we can print the median run time, highest 10%, lowest 10% etc
        //Maybe this can be an option we pass in? What degree of timing to use?
    //TODO: Make this work with multiple threads
public:
    TimingDetails() : m_indexClaimed(false), m_runCount(0), m_runTimesNanos(0.0) {}
    
    void claim(const std::string& _name) {
        if(m_indexClaimed && _name != m_name) {
            printAndThrowRuntimeError(_name + " is trying to claim a timing detail that is already claimed by " + m_name);
        }
        else {
            m_indexClaimed = true;
            m_name = _name;
        }
    }
    
    void clearClaim() {
        m_name = "";
        m_indexClaimed = false;
        m_runCount = 0;
        m_runTimesNanos = 0;
    }
    
    void stop() {
        m_end = std::chrono::high_resolution_clock::now();
        m_runTimesNanos += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_start).count());
        m_runCount++;
    }
    
    void start() {
        m_start = std::chrono::high_resolution_clock::now();
    }
        
    std::string getTimeString() const {
        return "Total Time: " + formatTimePretty(m_runTimesNanos) + "\tAverage Time: " + formatTimePretty(m_runTimesNanos / m_runCount);
    }
    
    std::string getFullString() const {
        return m_name + ":\tRun Count: " + std::to_string(m_runCount) + "\t" + getTimeString();
    }
    
    friend std::ostream& operator<<(std::ostream& stream, const TimingDetails& timingDetails) {
        stream << timingDetails.getFullString();
        return stream;
    }
    
private:
    std::string                                                 m_name;
    bool                                                        m_indexClaimed;
    size_t                                                      m_runCount;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
    double                                                      m_runTimesNanos;
};

class Timer
{
public:
    static Timer& getInstance()
    {
        static Timer    instance;
        return instance;
    }
    
    //Deleted copy functions
    Timer(Timer const&)           = delete;
    void operator=(Timer const&)  = delete;

private:
    Timer() {}
    
public:
    inline void startTimer(size_t index) {
#ifdef USE_TIMING
        if(!Timer::isEnabled()) {
            return;
        }
        m_timingDetails[index].start();
#endif
    }

    inline void stopTimer(size_t index) {
#ifdef USE_TIMING
        if(!Timer::isEnabled()) {
            return;
        }
        m_timingDetails[index].stop();
#endif
    }
    
    inline void registerTimer(size_t& index, std::string name) {
#ifdef USE_TIMING
        if(!Timer::isEnabled()) {
            return;
        }
        
        constexpr size_t UnusedIndexVal = -1;
        if(index == UnusedIndexVal) {
            index = m_index++;
            m_timingDetails.resize(m_index);
        }
                
        m_timingDetails[index].claim(name);
#endif
    }
    
    inline static void getTimingResultsAndClear() {
#ifdef USE_TIMING
        //Only print the results if timing is enabled
        if(m_enabled) {
#ifdef TESTING
            getInstance().printTimingResults();
#else
            getInstance().logResults();
#endif
        }
        getInstance().clearClaims();
#endif
    }
    
    static void enable() {
        m_enabled = true;
    }
    
    static void disable() {
        m_enabled = false;
    }
    
    static bool isEnabled() {
        return m_enabled;
    }
    
private:
    
    inline void printTimingResults() {
        printResults(std::cout);
    }
    
    void logResults() {
        std::ofstream file;
        file.open ("timing.txt");
        printResults(file);
        file.close();

    }
    
    void printResults(std::ostream& stream) { //TODO: Make this sorted by time spent.
        stream<<"\nTIMING RESULTS\n";
        for(size_t i = 0; i < m_timingDetails.size(); i++) {
            stream<<m_timingDetails[i]<<"\n";
        }
        stream<<"\n";
    }
    
    void clearClaims() {
        for(size_t i = 0; i < m_timingDetails.size(); i++) {
            m_timingDetails[i].clearClaim();
        }
    }
    
    static size_t                   m_index;
    static bool                     m_enabled;
    std::vector<TimingDetails>      m_timingDetails;
};
bool Timer::m_enabled = false;
size_t Timer::m_index = 0;

#endif /* Timer_h */
