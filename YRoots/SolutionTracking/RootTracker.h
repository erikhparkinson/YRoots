//
//  RootTracker.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef RootTracker_h
#define RootTracker_h

#include "utilities.h"
#include <complex>

struct FoundRoot {
    std::vector<double> root;
    SolveMethod         solveMethod;
    Interval            interval;
    double              conditionNumber;
};

class RootTracker {
public:
    RootTracker(size_t _numThreads) : m_numThreads(_numThreads)
    {
        m_foundRoots.resize(m_numThreads);
    }
    
    void storeRoot(size_t threadNum, std::vector<std::complex<double>>& _root, Interval& _interval, SolveMethod _howFound, double _conditionNumber, double _goodZerosTol) {
        //Check if the root is in the boundary
        for(size_t i = 0; i < _root.size(); i++) {
            if(std::abs(std::real(_root[i])) > 1 + _goodZerosTol) {
                return;
            }
            else if (std::abs(std::imag(_root[i])) > _goodZerosTol) {
                return;
            }
        }
        
        //Created the stored root.
        //TODO: Make this more efficient.
        m_foundRoots[threadNum].resize(m_foundRoots[threadNum].size() + 1);
        FoundRoot& thisRoot = m_foundRoots[threadNum][m_foundRoots[threadNum].size() - 1];
        
        
        //Store the root information
        thisRoot.conditionNumber = _conditionNumber;
        thisRoot.interval = _interval;
        thisRoot.solveMethod = _howFound;

        //Transform the root
        thisRoot.root.resize(_root.size());
        for(size_t i = 0; i < _root.size(); i++) {
            thisRoot.root[i] = ((_interval.upperBounds[i] - _interval.lowerBounds[i]) * std::real(_root[i]) + (_interval.upperBounds[i] + _interval.lowerBounds[i])) /2.0;
        }
    }
    
    void printResults() {
        std::cout<<"Roots found:\n";
        for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
            for(size_t rootNum = 0; rootNum < m_foundRoots[threadNum].size(); rootNum++) {
                for(size_t i = 0; i < m_foundRoots[threadNum][rootNum].root.size(); i++) {
                    std::cout<<m_foundRoots[threadNum][rootNum].root[i]<<"\t";
                }
            }
            std::cout<<"\n";
        }
    }

private:
    size_t                               m_numThreads;
    std::vector<std::vector<FoundRoot>>  m_foundRoots;
};


#endif /* RootTracker_h */
