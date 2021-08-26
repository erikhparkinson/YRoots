//
//  RootTracker.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef RootTracker_h
#define RootTracker_h

#include <complex>
#include <fstream>
#include <iomanip>
#include "Utilities/utilities.hpp"
#include "Utilities/ErrorTracker.hpp"
#include "Functions/Function.hpp"

struct FoundRoot {
    std::vector<double> root;
    std::vector<double> residuals;
    std::vector<double> evalErrorsAtRoot;
    SolveMethod         solveMethod;
    Interval            interval;
    double              conditionNumber;
    
    bool operator < (const FoundRoot &otherRoot)
    {
        for(size_t i = 0; i < root.size(); i++) {
            if(root[i] < otherRoot.root[i]) {
                return true;
            }
        }
        return false;
    }

    bool operator == (const FoundRoot &otherRoot)
    {
        for(size_t i = 0; i < root.size(); i++) {
            if(root[i] != otherRoot.root[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator != (const FoundRoot &otherRoot)
    {
        return !(*this == otherRoot);
    }

    friend std::ostream& operator<<(std::ostream& stream, const FoundRoot& foundRoots) {
        int precision = std::numeric_limits<double>::digits10 + 1;
        
        for(size_t i = 0; i < foundRoots.root.size(); i++) {
            stream<<std::setprecision(precision)<<foundRoots.root[i];
            if(i+1 < foundRoots.root.size()){
                stream<<"\t";
            }
        }
        return stream;
    }
};

class RootTracker {
public:
    RootTracker(size_t _numThreads, std::vector<std::vector<Function::SharedFunctionPtr>>& _functions, const GeneralParameters& _generalParameters) :
    m_numThreads(_numThreads),
    m_allFunctions(_functions),
    m_computeResiduals(_generalParameters.computeResiduals)
    {
        m_foundRoots.resize(m_numThreads);
        m_outputFileRoots = "roots.csv";
        m_outputFileResiduals = "residuals.csv";
    }
    
    bool storeRoot(size_t threadNum, std::vector<std::complex<double> >& _root, Interval& _interval, SolveMethod _howFound, double _conditionNumber, double _goodZerosTol) {
        //Check if the root is in the boundary
        for(size_t i = 0; i < _root.size(); i++) {
            if(std::abs(std::real(_root[i])) > 1 + _goodZerosTol) {
                return false;
            }
            else if (std::abs(std::imag(_root[i])) > _goodZerosTol) {
                return false;
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
        
        //Compute the resisuals of the root
        if(m_computeResiduals) {
            evaluateResiduals(thisRoot, m_allFunctions[threadNum]);
        }

        return true;
    }

    void storeRoot(size_t threadNum, std::vector<double>& _root, Interval& _interval, SolveMethod _howFound, double _conditionNumber) {
        //For a root that is already transformed and real, and we know is good.
        
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
            thisRoot.root[i] = _root[i];
        }
        
        //Compute the resisuals of the root
        if(m_computeResiduals) {
            evaluateResiduals(thisRoot, m_allFunctions[threadNum]);
        }
    }
    
    void evaluateResiduals(FoundRoot& _currRoot, std::vector<Function::SharedFunctionPtr>& _functions) {
        for(size_t i = 0; i < _functions.size(); i++) {
            ErrorTracker result = _functions[i]->evaluate<ErrorTracker>(_currRoot.root);
            _currRoot.residuals.push_back(result.value);
            _currRoot.evalErrorsAtRoot.push_back(result.error);
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
        
    void logResults() {
        std::ofstream file;
        int precision = std::numeric_limits<double>::digits10 + 1;
        file.open (m_outputFileRoots);
        for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
            for(size_t rootNum = 0; rootNum < m_foundRoots[threadNum].size(); rootNum++) {
                for(size_t i = 0; i < m_foundRoots[threadNum][rootNum].root.size(); i++) {
                    file<<std::setprecision(precision)<<m_foundRoots[threadNum][rootNum].root[i];
                    if(i + 1 < m_foundRoots[threadNum][rootNum].root.size()){
                        file<<",";
                    }
                }
                if(rootNum + 1 < m_foundRoots[threadNum].size()){
                    file<<"\n";
                }
            }
            if(threadNum + 1 < m_numThreads){
                file<<"\n";
            }
        }
        file.close();
        
        if(m_computeResiduals) {
            logResiduals();
        }
    }

    void logResiduals() {
        std::ofstream file;
        int precision = std::numeric_limits<double>::digits10 + 1;
        file.open (m_outputFileResiduals);
        for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
            for(size_t rootNum = 0; rootNum < m_foundRoots[threadNum].size(); rootNum++) {
                for(size_t i = 0; i < m_foundRoots[threadNum][rootNum].root.size(); i++) {
                    file<<std::setprecision(precision)<<m_foundRoots[threadNum][rootNum].residuals[i];
                    file<<",";
                    file<<std::setprecision(precision)<<m_foundRoots[threadNum][rootNum].evalErrorsAtRoot[i];
                    if(i + 1 < m_foundRoots[threadNum][rootNum].root.size()){
                        file<<",";
                    }
                }
                if(rootNum + 1 < m_foundRoots[threadNum].size()){
                    file<<"\n";
                }
            }
            if(threadNum + 1 < m_numThreads){
                file<<"\n";
            }
        }
        file.close();
    }

    std::vector<FoundRoot> getRoots() {
        //printResults();
        std::vector<FoundRoot> allRoots;
        for(size_t threadNum = 0; threadNum < m_numThreads; threadNum++) {
            for(size_t rootNum = 0; rootNum < m_foundRoots[threadNum].size(); rootNum++) {
                allRoots.push_back(m_foundRoots[threadNum][rootNum]);
            }
        }
        return allRoots;
    }

private:
    size_t                                  m_numThreads;
    std::vector<std::vector<Function::SharedFunctionPtr>>& m_allFunctions;
    const bool                              m_computeResiduals;
    std::vector<std::vector<FoundRoot> >    m_foundRoots;
    std::string                             m_outputFileRoots;
    std::string                             m_outputFileResiduals;
};


#endif /* RootTracker_h */
