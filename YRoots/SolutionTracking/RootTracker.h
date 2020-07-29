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
#include <tbb/concurrent_vector.h>
#include <complex>

//THIS CLASS MUST BE THREAD SAFE

struct FoundRoot {
    std::vector<double> root;
    SolveMethod         solveMethod;
    Interval            interval;
    double              conditionNumber;
};

class RootTracker {
public:
    RootTracker() {
        
    }
    
    void storeRoot(std::vector<std::complex<double>>& _root, Interval& _interval, SolveMethod _howFound, double _conditionNumber, double _goodZerosTol) {
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
        tbb::concurrent_vector<FoundRoot>::iterator foundIter = m_foundRoots.grow_by(1);

        //Store the root information
        foundIter->conditionNumber = _conditionNumber;
        foundIter->interval = _interval;
        foundIter->solveMethod = _howFound;

        //Transform the root
        foundIter->root.resize(_root.size());
        for(size_t i = 0; i < _root.size(); i++) {
            foundIter->root[i] = ((_interval.upperBounds[i] - _interval.lowerBounds[i]) * std::real(_root[i]) + (_interval.upperBounds[i] + _interval.lowerBounds[i])) /2.0;
        }
    }
    
    void printResults() {
        std::cout<<"Roots found:\n";
        for(size_t rootNum = 0; rootNum < m_foundRoots.size(); rootNum++) {
            for(size_t i = 0; i < m_foundRoots[rootNum].root.size(); i++) {
                std::cout<<m_foundRoots[rootNum].root[i]<<"\t";
            }
            std::cout<<"\n";
        }
    }

private:
    tbb::concurrent_vector<FoundRoot>       m_foundRoots;
};


#endif /* RootTracker_h */
