//
//  Subdivision.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/27/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef Subdivision_Solver_h
#define Subdivision_Solver_h

#include "utilities.h"
#include "IntervalData.h"
#include "ChebyshevApproximator.h"

struct SubdivisionParameters {
    double relApproxTol = 1.e-15;
    double absApproxTol = 1.e-12;
    double maxConditionNumber = 1e5;
    double goodZerosFactor = 100;
    double minGoodZerosTol = 1e-5;
    bool checkEvaluationError = true;
    size_t checkEvaluationFrequency = 1;
    size_t approximationDegree = 3;
    size_t targetDegree = 10;
    size_t maxLevel = 999;
    bool returnPotentials = false;
    std::string method = "svd";
    double targetTol = 1.e-15;
    bool useTargetTol = true;
};

template <Dimension D>
class SubdivisionSolver
{
public:
    SubdivisionSolver(const std::vector<std::unique_ptr<FunctionInterface>>& _functions);
    
    void solve(Interval _currentInterval, size_t currentLevel);
    
private:
    const std::vector<std::unique_ptr<FunctionInterface>>& m_functions;
    SubdivisionParameters m_subdivisionParameters;
    IntervalData        m_intervalData;
    std::vector<std::unique_ptr<ChebyshevApproximator<D>>> m_chebyshevApproximators;
    
    
    //TODO: Create Interval Data Class
    //TODO: Create objects to handle root tracking
    //TODO: Pass in the tolerances
    //TODO: Make a class for FFTs
    //
    
    
};

#include "SubdivisionSolver1D.ipp"
#include "SubdivisionSolverND.ipp"

#endif /* Subdivision_Solver_h */
