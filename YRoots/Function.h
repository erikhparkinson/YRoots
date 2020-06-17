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
    
    void evaluateGrid(const std::vector<std::vector<double>>& grid, double* results)
    {
        //Do nothing in the 0-dimensional case
        size_t dimension = grid.size();
        if(dimension == 0) {
            return;
        }
        
        //Set up the needed variables
        size_t evalSpot = 0;
        size_t numPoints = grid[0].size();
        std::vector<double> inputPoints(dimension, 0.0);
        std::vector<size_t> inputSpot(dimension);
        for(size_t i = 0; i < dimension; i++) {
            inputPoints[i] = grid[i][0];
        }
        
        //Iterate through all the combinations
        size_t spotToInc = 0;
        results[evalSpot++] = evaluate(inputPoints);
        while (spotToInc < dimension) {
            bool firstPass = true;
            while(++inputSpot[spotToInc] < numPoints) {
                inputPoints[spotToInc] = grid[spotToInc][inputSpot[spotToInc]];
                results[evalSpot++] = evaluate(inputPoints);
                if(firstPass && spotToInc != 0) {
                    spotToInc = 0;
                }
                firstPass = false;
            }
            inputSpot[spotToInc] = 0;
            inputPoints[spotToInc] = grid[spotToInc][0];
            spotToInc++;
        }
    }
    
    void evaluatePoints(const std::vector<std::vector<double>>& inputPoints, double* results)
    {
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
