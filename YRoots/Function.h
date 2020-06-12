//
//  Function.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/23/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef Function_h
#define Function_h

#include <math.h>
#include "utilities.h"

class FunctionInterface{
public:
    FunctionInterface(std::string _functionString, std::vector<std::string> _variableNames):
    m_functionString(_functionString),
    m_variableNames(_variableNames),
    m_dimension(m_variableNames.size())
    {

    }
    
    virtual ~FunctionInterface() = default;
    
    virtual double evaluate(const std::vector<double>& inputPoints) = 0;
    
    void evaluatePoints(const std::vector<std::vector<double>>& inputPoints, double* results) {
        //TODO: Make a better call to this. Have a specific evaluate grid function.
        for(size_t i = 0; i < inputPoints.size(); i++) {
            results[i] = evaluate(inputPoints[i]);
        }
    }
    
    std::string to_string() {
        return m_functionString;
    }
    
    size_t getDimension() {
        return m_dimension;
    }
    
private:
    std::string                                 m_functionString;
    std::vector<std::string>                    m_variableNames;
    size_t                                      m_dimension;
};

#endif /* Function_h */
