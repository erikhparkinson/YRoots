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
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include "Utilities/utilities.hpp"
#include "Utilities/Timer.hpp"
#include "Functions/Polynomial.hpp"

enum FunctionType {
    SIN, //Syntax: sin(x)
    COS, //Syntax: cos(x)
    TAN, //Syntax: tan(x)
    COSH, //Syntax: cosh(x)
    SINH, //Syntax: sinh(x)
    TANH, //Syntax: tanh(x)
    LOG, //Syntax: log(x)
    LOG10, //Syntax: log10(x)
    LOG2, //Syntax: log2(x)
    POWER, //Syntax: x^y
    SQRT, //Syntax: sqrt(x)
    EXP, //Syntax: exp(x)
    CONSTANT, //Syntax: 4, e, pi
    CHEBYSHEV, //Syntax: T5(...), -7T435(...)
    POWER_BASIS_MONOMIAL, //Syntax: x^2, -2x, 2x^2*y^7z
    CHEBYSHEV_BASIS_MONOMIAL, //Syntax: T5(x), -7T5(x)*T3(y)
    //Used in Evaluation Structure
    SUM, //A sum of subfunctions + or -
    PRODUCT, //A product of subfunctions, * or /
    POWER_BASIS_POLYNOMIAL, //Sums of POWER_BASIS_MONOMIAL, for faster evaluation.
    CHEBYSHEV_BASIS_POLYNOMIAL, //Sums of CHEBYSHEV_BASIS_MONOMIAL, for faster evaluation.
    VARIABLE, //A single monomial of one variable. Includes signed coefficient.
    UNKNOWN //Something went wrong
};

//TOOD: If it fails to parse make it call a class error message that prints the name of the function and detials.
//      That way we always can see the name of a function that it can't parse.

enum EvaluateGridType {
    BASE,
    SIMPLE,
    COMBINE
};

struct CHAR {
    static constexpr char LEFT_PAREN = '(';
    static constexpr char RIGHT_PAREN = ')';

    static constexpr char PLUS = '+';
    static constexpr char MINUS = '-';
    static constexpr char DIVIDE = '/';
    static constexpr char TIMES = '*';
    static constexpr char POWER = '^';

    static constexpr char COMMA = ',';
};

struct EvaluateGridInfo {
    bool precomputed = false;
    std::vector<std::vector<size_t> > childEvalIndexes;
    size_t childEvalSize;
    
    //For Copying over the top level of top functions that aren't full dimension
    std::vector<size_t> topResultMap;
};

#define SIGNCHECKSUM(isPositive, result, number) (isPositive ? result += number : result -= number)
#define SIGNCHECKPRODUCT(isPositive, result, number) (isPositive ? result *= number : result /= number)

class Function {
public:
    typedef std::shared_ptr<Function> SharedFunctionPtr;
    typedef std::unordered_map<std::string, SharedFunctionPtr> FunctionMap;

    //The _functionString is copied in the constructor because it get's modified in the function parsing
    Function(const std::string& _functionName, std::string _functionString, const std::vector<std::string>& _variableNames) :
    m_isTopFunction(false),
    m_functionLevel(0),
    m_functionType(FunctionType::UNKNOWN),
    m_value(0),
    m_varIndex(0),
    m_functionString(_functionString),
    m_functionName(_functionName),
    m_variableNames(_variableNames),
    m_polynomial(_variableNames.size()),
    m_rank(_variableNames.size())
    {
        //Parse the function
        functionParse(m_functionString);
    }
    
    void defineFunctionDimensions() {
        //Populate m_hasDimension
        m_hasDimension.resize(m_rank);
        switch(m_functionType) {
            case FunctionType::CONSTANT:
                break;
            case FunctionType::VARIABLE:
                m_hasDimension[m_varIndex] = true;
                break;
            case FunctionType::POWER_BASIS_MONOMIAL:
                m_monomial.prepEvaluation();
                m_hasDimension = m_monomial.getHasDimension();
                break;
            case FunctionType::CHEBYSHEV_BASIS_MONOMIAL: //TODO: Write this for Chebyshev Monomials
                break;
            case FunctionType::POWER_BASIS_POLYNOMIAL:
                m_polynomial.prepEvaluation();
                m_hasDimension = m_polynomial.getHasDimension();
                break;
            case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL: //TODO: Write this for Chebyshev Polynomials
                break;
            default: {
                //Everything that is a base function should be in the case statements above.
                //Everything else uses it's subfunctions
                for(size_t i = 0; i < m_subfunctions.size(); i++) {
                    const std::vector<bool>& subDimensions = m_subfunctions[i]->getHasDimension();
                    for(size_t j = 0; j < subDimensions.size(); j++) {
                        m_hasDimension[j] = m_hasDimension[j] || subDimensions[j];
                    }
                }
                break;
            }
        }
                
        //The number of true spots in m_hasDimension is m_numUsedDimensions
        m_numUsedDimensions = 0;
        for(size_t i = 0; i < m_hasDimension.size(); i++) {
            if(m_hasDimension[i]) {
                m_numUsedDimensions++;
            }
        }
        
        //Define the m_evaluateGridType
        switch(m_functionType) {
            case FunctionType::CONSTANT:
            case FunctionType::VARIABLE:
            case FunctionType::POWER_BASIS_MONOMIAL:
            case FunctionType::CHEBYSHEV_BASIS_MONOMIAL:
            case FunctionType::POWER_BASIS_POLYNOMIAL:
            case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL:
                m_evaluateGridType = EvaluateGridType::BASE;
                break;
            case FunctionType::SIN:
            case FunctionType::COS:
            case FunctionType::TAN:
            case FunctionType::SINH:
            case FunctionType::COSH:
            case FunctionType::TANH:
            case FunctionType::LOG:
            case FunctionType::LOG10:
            case FunctionType::LOG2:
            case FunctionType::SQRT:
            case FunctionType::EXP:
            case FunctionType::CHEBYSHEV:
                m_evaluateGridType = EvaluateGridType::SIMPLE;
                break;
            case FunctionType::POWER:
            case FunctionType::SUM:
            case FunctionType::PRODUCT:
                m_evaluateGridType = EvaluateGridType::COMBINE;
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid main! Fix Switch Statement!");
                break;
        }
    }
    
    void prepEvaluatGrid(size_t _gridSize, EvaluateGridInfo& _infoToPopulate) {
        _infoToPopulate.precomputed = true;
        
        //Get the childEvalSize
        _infoToPopulate.childEvalSize = power(_gridSize, m_numUsedDimensions);

        //Make sure the m_partialEvaluations is big enough
        if(m_partialEvaluations.size() <= _infoToPopulate.childEvalSize) {
            m_partialEvaluations.resize(_infoToPopulate.childEvalSize);
        }

        if(m_evaluateGridType == EvaluateGridType::COMBINE) {
            //We need to create the childEvalIndexes for the non simple cases
            _infoToPopulate.childEvalIndexes.resize(_infoToPopulate.childEvalSize);
            
            //Start on the first dimensions
            std::vector<size_t> tempDimCount(m_rank);
            size_t dimToInc = 0;
            while(!m_hasDimension[dimToInc] && dimToInc < m_rank) {
                dimToInc++;
            }
            
            size_t firstUsedDim = dimToInc;
            size_t currIndexSpot = 0;
            //Iterate through the combinations
            bool finishedLoop = false;
            while (!finishedLoop) {
                //Iterate in the current dimension
                while(tempDimCount[dimToInc] < _gridSize) {
                    //Get the points for this spot
                    for(size_t subfunctionNum = 0; subfunctionNum < m_subfunctions.size(); subfunctionNum++) {
                        const std::vector<bool>& subDimensions =  m_subfunctions[subfunctionNum]->getHasDimension();
                        size_t indexToUse = 0;
                        size_t dimMultiplier = 1;
                        for(size_t currDim = 0; currDim < subDimensions.size(); currDim++) {
                            if(subDimensions[currDim]) {
                                indexToUse += dimMultiplier * tempDimCount[currDim];
                                dimMultiplier *= _gridSize;
                            }
                        }
                        
                        _infoToPopulate.childEvalIndexes[currIndexSpot].push_back(indexToUse);
                    }
                    
                    //Continue the loop
                    currIndexSpot++;
                    tempDimCount[dimToInc]++;
                }
                
                //Increment the dimension
                do {
                    //Zero out the current spot and continue
                    tempDimCount[dimToInc++] = 0;
                    //Go the the next spot
                    while(!m_hasDimension[dimToInc] && dimToInc < m_rank) {
                        dimToInc++;
                    }
                    //Check if it overflows as well
                    if(dimToInc < m_rank) {
                        tempDimCount[dimToInc]++;
                    }
                    else {
                        finishedLoop = true;
                        break;
                    }
                }
                while(tempDimCount[dimToInc] == _gridSize);
                dimToInc = firstUsedDim;
            }
        }
        
        if(m_isTopFunction && m_numUsedDimensions != m_rank) {
            const size_t numResultPoints = power(_gridSize, m_rank);
            _infoToPopulate.topResultMap.resize(numResultPoints);
            
            //Start on the first dimension
            std::vector<size_t> tempDimCount(m_rank);
            size_t dimToInc = 0;
            
            size_t currIndexSpot = 0;
            //Iterate through the combinations
            bool finishedLoop = false;
            while (!finishedLoop) {
                //Iterate in the current dimension
                while(tempDimCount[dimToInc] < _gridSize) {
                    size_t indexToUse = 0;
                    size_t dimMultiplier = 1;
                    for(size_t currDim = 0; currDim < m_hasDimension.size(); currDim++) {
                        if(m_hasDimension[currDim]) {
                            indexToUse += dimMultiplier * tempDimCount[currDim];
                            dimMultiplier *= _gridSize;
                        }
                    }
                    _infoToPopulate.topResultMap[currIndexSpot] = indexToUse;

                    //Continue the loop
                    currIndexSpot++;
                    tempDimCount[dimToInc]++;
                }
                
                //Increment the dimension
                do {
                    //Zero out the current spot and continue
                    tempDimCount[dimToInc++] = 0;
                    //Go the the next spot
                    if(dimToInc < m_rank) {
                        tempDimCount[dimToInc]++;
                    }
                    else {
                        finishedLoop = true;
                        break;
                    }
                }
                while(tempDimCount[dimToInc] == _gridSize); //Check if the next spot overflows
                dimToInc = 0;
            }
        }
    }
    
    void evaluateGridBase(const std::vector<std::vector<double> >& _grid, std::vector<double>& _results) {
        switch(m_functionType) {
            case FunctionType::CONSTANT:
                _results[0] =  m_value;
                break;
            case FunctionType::VARIABLE:
                for(size_t i = 0; i < _grid[0].size(); i++) {
                    _results[i] =  m_value * _grid[m_varIndex][i];
                }
                break;
            case FunctionType::POWER_BASIS_MONOMIAL:
                m_monomial.evaluateGrid(_grid, _results);
                break;
            //case FunctionType::CHEBYSHEV_BASIS_MONOMIAL:  //TODO: Write the specialized base call for CHEBYSHEV_BASIS_MONOMIAL
            //    break;
            case FunctionType::POWER_BASIS_POLYNOMIAL:
                m_polynomial.evaluateGrid(_grid, _results);
                break;
            //case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL:  //TODO: Write the specialized base call for CHEBYSHEV_BASIS_POLYNOMIAL
            //    break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid base! Fix Switch Statement!");
                break;
        }
    }

    void evaluateGridSimple(size_t _gridSize, std::vector<double>& _results) {
        //Get the childs evaluations and resize the result vector if needed
        const std::vector<double>& childEvals = m_subfunctions[0]->getPartialEvals();
        size_t numEvals = m_evaluateGridInfo[_gridSize].childEvalSize;

        //Run the functions over each evaluation
        switch(m_functionType) {
            case FunctionType::SIN:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * sin(childEvals[i]);
                }
                break;
            case FunctionType::COS:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * cos(childEvals[i]);
                }
                break;
            case FunctionType::TAN:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * tan(childEvals[i]);
                }
                break;
            case FunctionType::SINH:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * sinh(childEvals[i]);
                }
                break;
            case FunctionType::COSH:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * cosh(childEvals[i]);
                }
                break;
            case FunctionType::TANH:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * tanh(childEvals[i]);
                }
                break;
            case FunctionType::LOG:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * log(childEvals[i]);
                }
                break;
            case FunctionType::LOG10:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * log10(childEvals[i]);
                }
                break;
            case FunctionType::LOG2:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * log2(childEvals[i]);
                }
                break;
            case FunctionType::SQRT:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * sqrt(childEvals[i]);
                }
                break;
            case FunctionType::EXP:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * exp(childEvals[i]);
                }
                break;
            case FunctionType::CHEBYSHEV:
                for(size_t i = 0; i < numEvals; i++) {
                    //TODO: This is inefficient, make a function that does all of them at once.
                    _results[i] =  m_value * chebPower(childEvals[i], m_varIndex);
                }
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid main! Fix Switch Statement!");
                break;
        }
    }

    void evaluateGridCombine(size_t _gridSize, std::vector<double>& _results) {
        std::vector<std::vector<size_t> >& currSpots = m_evaluateGridInfo[_gridSize].childEvalIndexes;
        switch(m_functionType) {
            case FunctionType::POWER: {
                const std::vector<double>& childEvals1 = m_subfunctions[0]->getPartialEvals();
                const std::vector<double>& childEvals2 = m_subfunctions[1]->getPartialEvals();
                for(size_t i = 0; i < currSpots.size(); i++) {
                    _results[i] =  m_value * pow(childEvals1[currSpots[i][0]], childEvals2[currSpots[i][1]]);
                }
                break;
            }
            case FunctionType::SUM:{
                double result;
                for(size_t i = 0; i < currSpots.size(); i++) {
                    result = m_value;
                    for(size_t j = 0; j < currSpots[i].size(); j++) {
                        SIGNCHECKSUM(m_operatorSigns[j], result, m_subfunctions[j]->getPartialEvals()[currSpots[i][j]]);
                    }
                    _results[i] = result;
                }
                break;
            }
            case FunctionType::PRODUCT:{
                double result;
                for(size_t i = 0; i < currSpots.size(); i++) {
                    result = m_value;
                    for(size_t j = 0; j < currSpots[i].size(); j++) {
                        SIGNCHECKPRODUCT(m_operatorSigns[j], result, m_subfunctions[j]->getPartialEvals()[currSpots[i][j]]);
                    }
                    _results[i] = result;
                }
                break;
            }
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid main! Fix Switch Statement!");
                break;
        }
    }
    
    void evaluateGridMain(const std::vector<std::vector<double> >& _grid) {
        //Prep the evaluate grid info if it hasn't already been done
        size_t gridSize = _grid[0].size();
        if(unlikely(m_evaluateGridInfo.size() <= gridSize)) {
            m_evaluateGridInfo.resize(gridSize + 1);
        }
        if(unlikely(!m_evaluateGridInfo[gridSize].precomputed)) {
            prepEvaluatGrid(gridSize, m_evaluateGridInfo[gridSize]);
        }

        switch(m_evaluateGridType) {
            case EvaluateGridType::BASE:
                evaluateGridBase(_grid, m_partialEvaluations);
                break;
            case EvaluateGridType::SIMPLE:
                evaluateGridSimple(gridSize, m_partialEvaluations);
                break;
            case EvaluateGridType::COMBINE:
                evaluateGridCombine(gridSize, m_partialEvaluations);
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid main! Fix Switch Statement!");
                break;
        }
    }

    void evaluateGridMain(const std::vector<std::vector<double> >& _grid, std::vector<double>& _results) {
        assert(m_isTopFunction); //This should only be called by the top function.
        //Prep the evaluate grid info if it hasn't already been done
        size_t gridSize = _grid[0].size();
        if(unlikely(m_evaluateGridInfo.size() <= gridSize)) {
            m_evaluateGridInfo.resize(gridSize + 1);
        }
        if(unlikely(!m_evaluateGridInfo[gridSize].precomputed)) {
            prepEvaluatGrid(gridSize, m_evaluateGridInfo[gridSize]);
        }
        
        //If we are not using all the dimensions, run the evaluations with m_partialEvaluations, then copy it into _results after
        const bool usesAllDims = m_numUsedDimensions == m_rank;
        std::vector<double>& resultSpot = usesAllDims ? _results : m_partialEvaluations;
        
        switch(m_evaluateGridType) {
            case EvaluateGridType::BASE:
                evaluateGridBase(_grid, resultSpot);
                break;
            case EvaluateGridType::SIMPLE:
                evaluateGridSimple(gridSize, resultSpot);
                break;
            case EvaluateGridType::COMBINE:
                evaluateGridCombine(gridSize, resultSpot);
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid main! Fix Switch Statement!");
                break;
        }
        
        if(!usesAllDims) { //Copy from m_partialEvaluations into _results.
            //For every dim that isn't used,
            for(size_t i = 0; i < m_evaluateGridInfo[gridSize].topResultMap.size(); i++) {
                _results[i] = m_partialEvaluations[m_evaluateGridInfo[gridSize].topResultMap[i]];
            }
        }
    }
    
    void evaluateGrid(const std::vector<std::vector<double> >& _grid, std::vector<double>& _results) {
        //Evaluates the grid so that the evaluation of (grid[0][i], grid[1][j]) is in _results[j+grid.size() + i]
        
        //Create the function tree to evaluate.
        if(unlikely(!m_isTopFunction)) {
            m_isTopFunction = true;
            getFunctionLevels(m_allFunctionLevels);
        }
        
        //Call everything level by level
        for(size_t level = 0; level < m_allFunctionLevels.size(); level++) {
            for(auto&& func: m_allFunctionLevels[level]) {
                func->evaluateGridMain(_grid);
            }
        }
        
        //Now call this top function.
        evaluateGridMain(_grid, _results);
        return;

        //This is the old code, it can still be run with the assert statement to  gives the same results.
        
        //Do nothing in the 0-dimensional case
        size_t dimension = _grid.size();
        if(dimension == 0) {
            return;
        }
        
        //Set up the needed variables
        size_t evalSpot = 0;
        size_t numPoints = _grid[0].size();
        std::vector<double> inputPoints(dimension, 0.0);
        std::vector<size_t> inputSpot(dimension);
        for(size_t i = 0; i < dimension; i++) {
            inputPoints[i] = _grid[i][0];
        }
        
        //Iterate through all the combinations
        size_t spotToInc = 0;
        //assert(_results[evalSpot++] == evaluate<double>(inputPoints));
        _results[evalSpot++] = evaluate<double>(inputPoints);
        while (spotToInc < dimension) {
            while(++inputSpot[spotToInc] < numPoints) {
                inputPoints[spotToInc] = _grid[spotToInc][inputSpot[spotToInc]];
                double diff = std::abs(_results[evalSpot] - evaluate<double>(inputPoints));
                double relativeDiff = diff / std::abs(_results[evalSpot]);
                //assert(relativeDiff < 1e-5); evalSpot++;
                _results[evalSpot++] = evaluate<double>(inputPoints);
                if(spotToInc != 0) {
                    spotToInc = 0;
                }
            }
            inputSpot[spotToInc] = 0;
            inputPoints[spotToInc] = _grid[spotToInc][0];
            spotToInc++;
        }
    }

    template<typename ReturnType>
    ReturnType evaluate(const std::vector<double>& _inputPoints) {
        switch(m_functionType) {
            case FunctionType::SIN:
                return static_cast<ReturnType>(m_value) * sin(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::COS:
                return static_cast<ReturnType>(m_value) * cos(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::TAN:
                return static_cast<ReturnType>(m_value) * tan(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::SINH:
                return static_cast<ReturnType>(m_value) * sinh(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::COSH:
                return static_cast<ReturnType>(m_value) * cosh(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::TANH:
                return static_cast<ReturnType>(m_value) * tanh(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::LOG:
                return static_cast<ReturnType>(m_value) * log(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::LOG10:
                return static_cast<ReturnType>(m_value) * log10(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::LOG2:
                return static_cast<ReturnType>(m_value) * log2(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::POWER:
                return static_cast<ReturnType>(m_value) * pow(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints), m_subfunctions[1]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::SQRT:
                return static_cast<ReturnType>(m_value) * sqrt(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::EXP:
                return static_cast<ReturnType>(m_value) * exp(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints));
            case FunctionType::CONSTANT:
                return static_cast<ReturnType>(m_value);
            case FunctionType::POWER_BASIS_POLYNOMIAL:
                return m_polynomial.evaluate<ReturnType>(_inputPoints);
            //case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL:
            //    break;
            case FunctionType::SUM:
                return sumEval<ReturnType>(_inputPoints);
            case FunctionType::PRODUCT:
                return productEval<ReturnType>(_inputPoints);
            case FunctionType::POWER_BASIS_MONOMIAL:
                return m_monomial.evaluate<ReturnType>(_inputPoints);
            //case FunctionType::CHEBYSHEV_BASIS_MONOMIAL:
            //    return chebyshevMonomialEval<ReturnType>(_inputPoints);
            case FunctionType::CHEBYSHEV:
                return static_cast<ReturnType>(m_value) * chebPower(m_subfunctions[0]->evaluate<ReturnType>(_inputPoints), m_varIndex);
            case FunctionType::VARIABLE:
                return static_cast<ReturnType>(m_value) * ((ReturnType)_inputPoints[m_varIndex]);
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate! Fix Switch Statement!");
                break;
        }
        return 0.0;
    }
        
private:
    void functionParse(std::string _functionString) { //Don't pass by reference, as we change it.
        //GRAMMER:
        // ( and ) make subfunctions
        // +,-,*,/,^ are math symbols
        //numbers, e, pi are constants
        //Everything else must be a FunctionType, string variable name, or subfunction name
                
        //Remove excess parenthesis surrounding the entire function.
        removeExtraParenthesis(_functionString);
        
        //To start, split on +-.
        std::vector<std::string> subparts;
        splitSum(_functionString, subparts, m_operatorSigns);
        
        if(subparts.size() == 0) { //This is an empty string! Something went wrong
            printAndThrowRuntimeError("Error Parsing Functions! Found Empty string in recursive step!");
        }
        //If m_operatorSigns[0] is false, then the function is of the form -<expression>. This must be parsed as a sum type.
        else if(subparts.size() == 1 && m_operatorSigns[0]) {
            //Reparse to split on */
            subparts.clear();
            m_operatorSigns.clear();
            splitProduct(_functionString, subparts, m_operatorSigns);
            
            if(subparts.size() == 1) {
                //Reparse to splt on ^
                subparts.clear();
                m_operatorSigns.clear();
                splitPower(_functionString, subparts);

                if(subparts.size() == 1) { //This is a more complex type.
                    parseComplexType(_functionString);
                }
                else { //Make this a POWER type and split it up
                    parsePower(subparts);
                }
            }
            else { //Make this a PRODUCT type and split it up
                parseProduct(subparts);
            }
        }
        else { //Make this a SUM type and split it up
            parseSum(subparts);
        }
        
        //Make simplifications that will make calculations easier in the future.
        simplifyExpressions();
    }
        
    void simplifyExpressions() {
        //TODO: To Update everything to shared ptrs, it is important to never change a function besdides the
        //one you are currently in. Otherwise it will be changing other copies of the function that are unrelated.
        
        //TODO: I think the only place that is an issue is with sums and products. So if a function appears the same
        //twice in a sum I should remove both and make a new function that is 2*<that function> (parens needed?)
        //If a function appears twice in a product then I should remove both and make a new function that is <that function>^2
        
        //TODO: If I have a monomial raised to a power I should be able to simplify that.
        
        
        //Pull any constants out of a sum
        if(m_functionType == FunctionType::SUM) {
            //Pull any constants out of a sum
            size_t funcIndex = 0;
            while(funcIndex < m_subfunctions.size()) {
                if(m_subfunctions[funcIndex]->getFunctionType() == FunctionType::CONSTANT) {
                    SIGNCHECKSUM(m_operatorSigns[funcIndex], m_value, m_subfunctions[funcIndex]->getValue());
                    m_subfunctions.erase(m_subfunctions.begin() + funcIndex);
                    m_operatorSigns.erase(m_operatorSigns.begin() + funcIndex);
                }
                else {
                    funcIndex++;
                }
            }
            //If the size of the sum is 0 it is now a constant
            if(m_subfunctions.size() == 0) {
                m_functionType = FunctionType::CONSTANT;
            }
        }
        
        //Pull any constants out of a product. Check if we multiply or divide!
        if(m_functionType == FunctionType::PRODUCT) {
            //Pull any constants out of a product
            size_t funcIndex = 0;
            while(funcIndex < m_subfunctions.size()) {
                if(m_subfunctions[funcIndex]->getFunctionType() == FunctionType::CONSTANT) {
                    SIGNCHECKPRODUCT(m_operatorSigns[funcIndex], m_value, m_subfunctions[funcIndex]->getValue());
                    m_subfunctions.erase(m_subfunctions.begin() + funcIndex);
                    m_operatorSigns.erase(m_operatorSigns.begin() + funcIndex);
                }
                else {
                    funcIndex++;
                }
            }
            //If the size of the product is 0 it is now a constant
            if(m_subfunctions.size() == 0) {
                m_functionType = FunctionType::CONSTANT;
            }
        }
        
        //Check for a constant LOG
        if(m_functionType == FunctionType::LOG) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= log(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a constant LOG10
        if(m_functionType == FunctionType::LOG10) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= log10(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a constant LOG2
        if(m_functionType == FunctionType::LOG2) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= log2(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant SIN
        if(m_functionType == FunctionType::SIN) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= sin(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant COS
        if(m_functionType == FunctionType::COS) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= cos(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant TAN
        if(m_functionType == FunctionType::TAN) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= tan(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant SINH
        if(m_functionType == FunctionType::SINH) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= sinh(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant COSH
        if(m_functionType == FunctionType::COSH) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= cosh(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant TANH
        if(m_functionType == FunctionType::TANH) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= tanh(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant SQRT
        if(m_functionType == FunctionType::SQRT) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= sqrt(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant EXP
        if(m_functionType == FunctionType::EXP) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= exp(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }

        //Check for a contant POWER
        if(m_functionType == FunctionType::POWER) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT && m_subfunctions[1]->getFunctionType() == FunctionType::CONSTANT) {
                m_value *= pow(m_subfunctions[0]->getValue(), m_subfunctions[1]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
            }
        }
        
        //Check powers with ones and zeros
        if(m_functionType == FunctionType::POWER) {
            //Check if it is something to the POWER of 0
            if(m_subfunctions[1]->getFunctionType() == FunctionType::CONSTANT && m_subfunctions[1]->getValue() == 0.0) {
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
                m_value = 1;
            }
            //Check if it is something to the POWER of 1
            else if(m_subfunctions[1]->getFunctionType() == FunctionType::CONSTANT && m_subfunctions[1]->getValue() == 1.0) {
                //TODO: Find a way to do this without copying. Maybe just reparse the subfunction string?
                
                //std::shared_ptr tempFunctionPtr(m_subfunctions[0]);
                //copyFunction(*tempFunctionPtr.get());
            }
            //Check it is 0 to the POWER of something
            else if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT && m_subfunctions[0]->getValue() == 0.0) {
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
                m_value = 0;
            }
            //Check it is 1 to the POWER of something
            else if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT && m_subfunctions[0]->getValue() == 1.0) {
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::CONSTANT;
                m_value = 1;
            }
        }

        //TODO: Think about how to do monomials and polynomials
        //return; //TODO: Remove this when we want to run polynomials
        
        //Check if a POWER is a POWER_BASIS_MONOMIAL
        if(m_functionType == FunctionType::POWER) {
            if(m_subfunctions[0]->getFunctionType() == FunctionType::VARIABLE && m_subfunctions[1]->getFunctionType() == FunctionType::CONSTANT) {
                size_t power = static_cast<size_t>(m_subfunctions[1]->getValue());
                m_monomial.clear(m_rank);
                if(power == m_subfunctions[1]->getValue()) {
                    m_monomial.spot[m_subfunctions[0]->getVarIndex()] = power;
                    m_monomial.coeff = m_subfunctions[0]->getValue();
                    m_subfunctions.erase(m_subfunctions.begin());
                    m_subfunctions.erase(m_subfunctions.begin());
                    m_functionType = FunctionType::POWER_BASIS_MONOMIAL;
                }
            }
        }
        
        //Check if a PRODUCT is composed of only POWER_BASIS_MONOMIALs, VARIABLEs and CONSTANTs. Then it should be just one big POWER_BASIS_MONOMIAL.
        if(m_functionType == FunctionType::PRODUCT) {
            //Check if everything is POWER_BASIS_MONOMIAL or VARIABLE
            bool isMonomial = true;
            bool isPolynomial = true;
            int polySpot = -1;
            for(size_t funcNum = 0; funcNum < m_subfunctions.size(); funcNum++) {
                const FunctionType subType = m_subfunctions[funcNum]->getFunctionType();
                const bool monStuff = subType == FunctionType::POWER_BASIS_MONOMIAL || subType == FunctionType::VARIABLE;
                const bool polyStuff = subType == FunctionType::POWER_BASIS_POLYNOMIAL || monStuff;
                isMonomial &= monStuff && m_operatorSigns[funcNum]; //Make sure it's not a division
                isPolynomial &= polyStuff && m_operatorSigns[funcNum]; //Make sure it's not a division
                if(subType == FunctionType::POWER_BASIS_POLYNOMIAL) {
                    if(polySpot == -1) { //Allow one polynomial
                        polySpot = funcNum;
                    }
                    else {
                        isPolynomial = false;
                    }
                }
            }
            //Create the new POWER_BASIS_MONOMIAL
            if(isMonomial) {
                m_monomial.clear(m_rank);
                m_monomial.coeff = m_value;
                for(size_t funcNum = 0; funcNum < m_subfunctions.size(); funcNum++) {
                    switch(m_subfunctions[funcNum]->getFunctionType()) {
                    case FunctionType::POWER_BASIS_MONOMIAL: {
                        //Update the coefficient.
                        m_monomial *= m_subfunctions[funcNum]->getMonomial();
                        break;
                    }
                    case FunctionType::VARIABLE: {
                        //Update the coefficient.
                        m_monomial.coeff *= m_subfunctions[funcNum]->getValue();
                        //Increment the index at the correct spot by 1
                        m_monomial.spot[m_subfunctions[funcNum]->getVarIndex()] += 1;
                        break;
                    }
                    case FunctionType::CONSTANT: {
                        //Update the coefficient.
                        m_monomial.coeff *= m_subfunctions[funcNum]->getValue();
                        break;
                    }
                    default: {
                        printAndThrowRuntimeError("Error converting product to Power Basis Monomial!");
                        break;
                    }
                    }
                }
                m_functionType = FunctionType::POWER_BASIS_MONOMIAL;
                m_subfunctions.clear();
            }
            else if(isPolynomial) {
                m_polynomial = m_subfunctions[polySpot]->getPolynomial();
                //Multiply by the constant term
                if(m_value != 1.0) {
                    Monomial toAdd;
                    toAdd.spot.resize(m_rank, 0);
                    toAdd.coeff = m_value;
                    m_polynomial.multiplyMonomial(toAdd);
                }
                //Multiply by everything else
                for(size_t funcNum = 0; funcNum < m_subfunctions.size(); funcNum++) {
                    if(funcNum == polySpot) {
                        continue;
                    }
                    switch(m_subfunctions[funcNum]->getFunctionType()) {
                    case FunctionType::POWER_BASIS_MONOMIAL: {
                        //Update the coefficient.
                        m_polynomial.multiplyMonomial(m_subfunctions[funcNum]->getMonomial());
                        break;
                    }
                    case FunctionType::VARIABLE: {
                        //Multiply this as a monomial
                        Monomial toAdd;
                        toAdd.spot.resize(m_rank, 0);
                        toAdd.spot[m_subfunctions[funcNum]->getVarIndex()] = 1;
                        toAdd.coeff = m_subfunctions[funcNum]->getValue();
                        m_polynomial.multiplyMonomial(toAdd);
                        break;
                    }
                    case FunctionType::CONSTANT: {
                        //Multiply this as a monomial
                        Monomial toAdd;
                        toAdd.spot.resize(m_rank, 0);
                        toAdd.coeff = m_subfunctions[funcNum]->getValue();
                        m_polynomial.multiplyMonomial(toAdd);
                        break;
                    }
                    default: {
                        printAndThrowRuntimeError("Error converting product to Power Basis Monomial!");
                        break;
                    }
                    }
                }
                m_functionType = FunctionType::POWER_BASIS_POLYNOMIAL;
                m_subfunctions.clear();
            }
        }
        
        //Check if a sum is composed of only POWER_BASIS_POLYNOMIALs, POWER_BASIS_MONOMIALs and VARIABLEs. The it is a POWER_BASIS_POLYNOMIAL.
        if(m_functionType == FunctionType::SUM) {
            //Check if everything is POWER_BASIS_MONOMIAL or VARIABLE
            bool isPolynomial = true;
            for(size_t funcNum = 0; funcNum < m_subfunctions.size(); funcNum++) {
                isPolynomial &= (m_subfunctions[funcNum]->getFunctionType() == FunctionType::POWER_BASIS_MONOMIAL || m_subfunctions[funcNum]->getFunctionType() == FunctionType::VARIABLE || m_subfunctions[funcNum]->getFunctionType() == FunctionType::POWER_BASIS_POLYNOMIAL);
            }
            if(isPolynomial) {
                //Add the monomials to the polynomial
                m_polynomial.clear();
                //Add the constant term.
                Monomial toAdd;
                if(m_value != 0.0) {
                    toAdd.spot.resize(m_rank, 0);
                    toAdd.coeff = m_value;
                    m_polynomial.addMonomial(toAdd);
                }
                //Add the monomials
                for(size_t funcNum = 0; funcNum < m_subfunctions.size(); funcNum++) {
                    switch(m_subfunctions[funcNum]->getFunctionType()) {
                    case FunctionType::POWER_BASIS_POLYNOMIAL: {
                        //Add all the monomials.
                        if(!m_operatorSigns[funcNum]) { //Flip the signs on the coeffs
                            for(const Monomial& m : m_subfunctions[funcNum]->getPolynomial().getMonomials()) {
                                toAdd = m;
                                toAdd.coeff *= -1;
                                m_polynomial.addMonomial(toAdd);
                            }
                        }
                        else {
                            m_polynomial.addMonomials(m_subfunctions[funcNum]->getPolynomial().getMonomials());
                        }
                        break;
                    }
                    case FunctionType::POWER_BASIS_MONOMIAL: {
                        //Add the monomial.
                        toAdd = m_subfunctions[funcNum]->getMonomial();
                        if(!m_operatorSigns[funcNum]) {
                            toAdd.coeff *= -1;
                        }
                        m_polynomial.addMonomial(toAdd);
                        break;
                    }
                    case FunctionType::VARIABLE: {
                        //Add it as a monomial.
                        toAdd.clear(m_rank);
                        toAdd.coeff = m_subfunctions[funcNum]->getValue() * (m_operatorSigns[funcNum] ? 1.0 : -1.0); //Get the sign based on the +-
                        toAdd.spot[m_subfunctions[funcNum]->getVarIndex()] = 1;
                        m_polynomial.addMonomial(toAdd);
                        break;
                    }
                    default: {
                        printAndThrowRuntimeError("Error converting product to Power Basis Monomial!");
                        break;
                    }
                    }
                }
                m_functionType = FunctionType::POWER_BASIS_POLYNOMIAL;
                m_subfunctions.clear();
            }
        }
        
        if(m_functionType == FunctionType::POWER_BASIS_MONOMIAL) {
            m_monomial.prepEvaluation();
        }
        
        //TODO: A product with 1 POWER_BASIS_POLYNOMIAL and POWER_BASIS_MONOMIALs and VARIABLEs
        
        
        //Check if a SUM has multiple POWER_BASIS_MONOMIAL, VARIABLE, and CONSTANT in it
        //TODO: Decide then this should be turned to a POWER_BASIS_POLYNOMIAL
        
        //The advantage of a POLYNOMIAL is that it does about 2 operations per monomial it holds
        //A sum of polynomials has to do the log reduction which is much slower, plus the function
        //Call and switch statement. Overall maybe 10 operations or so? But also probably slower ones.
        //I'd guess the change boundary is when it is between 1-10% full. Try something and then run timings
                
        //Check if a SUM has multiple CHEBYSHEV_BASIS_MONOMIAL and CONSTANT in it
        //TODO: Decide then this should be turned to a CHEBYSHEV_BASIS_POLYNOMIAL

        
        
        //TODO: Turn Sums of Monomials into Polynomials if that would make it more efficient
        //TODO: Check Sum Types to see if they should be polynomials!
        //Grab every element of the sum that is a monomial or chebyshev monomial and decide if it is worth pulling out.
        //Note: For a non-dense system it won't be. EX 1 + x^100 + y^100 is much better just being computed as three monomials then
        //Using Horners Method on it because of all of the 0's.
        //After pulling out the polys, if only normal or cheb poly exists, make that the type
        //Otherwise, it is sum type. Make the normal and cheb polys part of the sum
        //TODO: In that case make a seperate constructor that takes in the pre-parsed monomials
        
        
        //TODO: At the end go through and see if any Functions are equal to each other, if so make them the same shared ptr.
        //Make a Tree, constants and variables are level 0 and everything else is level 1 + max of sublevels.
        //When evaluating go through and evaluate each function in level 0 then 1, etc.
        //For evaluate grid pass in the size of the max vector needed. Then have each function store the values in a grid.
        
        //TODO: There should be no constants left at the end of this function! That should all be assimilated into the others.
        //Unless the whole function is one big constant I guess.
    }
        
    void removeExtraParenthesis(std::string& _functionString) {
        //Loop that removes parenthesis
        while(_functionString.length() > 0 && _functionString[0] == CHAR::LEFT_PAREN) {
            size_t parenthesisCount = 1;
            bool removeLayer = false;
            //Loop through the string looking for a space enclosed by parenthesis
            for (std::string::size_type stringSpot = 1; stringSpot < _functionString.length(); stringSpot++) {
                char currChar = _functionString[stringSpot];
                //Update the parenthesisCount
                if(currChar == CHAR::RIGHT_PAREN) {
                    parenthesisCount--;
                }
                else if(currChar == CHAR::LEFT_PAREN) {
                    parenthesisCount++;
                }
                //If this closes a parenthesis we are done
                if(parenthesisCount == 0) {
                    //Remove the parenthesis if it covers the whole thing
                    removeLayer = (stringSpot + 1 == _functionString.length());
                    break;
                }
            }
            //Remove a layer of parenthesis or stop.
            if(removeLayer) {
                _functionString = _functionString.substr(1, _functionString.length() - 2);
            }
            else {
                break;
            }
        }
    }
    
    std::unordered_set<size_t> getScientificNotationMinusLocs(const std::string& _functionString) {
        std::unordered_set<size_t> spots;
        for(size_t spot = 0; spot < _functionString.length(); spot++) {
            if(_functionString[spot] == CHAR::MINUS) {
                //It is scientific notation if it is <isNumericDigit><.>e-<.><isNumericDigit> with optional <.>
                const bool eBefore = spot >= 1 && (_functionString[spot-1] == 'e' || _functionString[spot-1] == 'E');
                const bool numericBefore = spot >= 2 && isNumericDigit(_functionString[spot-2]);
                const bool numericBefore2 = spot >= 3 && isNumericDigit(_functionString[spot-3]) && _functionString[spot-2] == '.';
                const bool numericAfter = _functionString.length()-spot-1 >= 1 && isNumericDigit(_functionString[spot+1]);
                const bool numericAfter2 = _functionString.length()-spot-1 >= 2 && _functionString[spot+1] == '.' && isNumericDigit(_functionString[spot+2]);
                if(eBefore && (numericBefore || numericBefore2) && (numericAfter || numericAfter2)) {
                    spots.insert(spot);
                }
            }
        }
        return spots;
    }
    
    void splitSum(const std::string& _functionString, std::vector<std::string>& _subparts, std::vector<bool>& _isSum) {
        std::unordered_set<size_t> scientificNotationMinusLocs = getScientificNotationMinusLocs(_functionString);
        
        std::string currentSubstring = "";
        bool currentIsSum = true;
        size_t parenthesisCount = 0;
        for (std::string::size_type stringSpot = 0; stringSpot < _functionString.size(); stringSpot++) {
            char currChar = _functionString[stringSpot];
            //Update the parenthesisCount
            if(currChar == CHAR::RIGHT_PAREN) {
                if(parenthesisCount == 0) {
                    printAndThrowRuntimeError("Mismatched Parenthesis in Function!");
                }
                parenthesisCount--;
                currentSubstring += currChar;
            }
            else if(currChar == CHAR::LEFT_PAREN) {
                parenthesisCount++;
                currentSubstring += currChar;
            }
            else if(parenthesisCount > 0) {
                //Do nothing, we are inside parenthesis
                currentSubstring += currChar;
            }
            else if(currChar == CHAR::MINUS && scientificNotationMinusLocs.find(stringSpot) != scientificNotationMinusLocs.end()) {
                //Do nothing, this MINUS is part of scientific notation.
                currentSubstring += currChar;
            }
            //Store the current substring if it exists
            else if(currChar == CHAR::PLUS || currChar == CHAR::MINUS) {
                if(currentSubstring.length() > 0) {
                    _subparts.push_back(currentSubstring);
                    _isSum.push_back(currentIsSum);
                }
                currentSubstring = "";
                currentIsSum = currChar == CHAR::PLUS;
            }
            else {
                currentSubstring += currChar;
            }
        }
        
        //Add Anything Remaining
        _subparts.push_back(currentSubstring);
        _isSum.push_back(currentIsSum);
    }

    void splitProduct(const std::string& _functionString, std::vector<std::string>& _subparts, std::vector<bool>& _isMultiply) {
        std::string currentSubstring = "";
        bool currentIsMultiply = true;
        size_t parenthesisCount = 0;
        for (std::string::size_type stringSpot = 0; stringSpot < _functionString.size(); stringSpot++) {
            char currChar = _functionString[stringSpot];
            //Update the parenthesisCount
            if(currChar == CHAR::RIGHT_PAREN) {
                if(parenthesisCount == 0) {
                    printAndThrowRuntimeError("Mismatched Parenthesis in Function!");
                }
                parenthesisCount--;
                currentSubstring += currChar;
            }
            else if(currChar == CHAR::LEFT_PAREN) {
                parenthesisCount++;
                currentSubstring += currChar;
            }
            else if(parenthesisCount > 0) {
                //Do nothing, we are inside parenthesis
                currentSubstring += currChar;

            }
            else if(currChar == CHAR::TIMES && stringSpot + 1 < _functionString.size() && _functionString[stringSpot+1] == CHAR::TIMES) {
                //Do nothing, this is a ** which is a ^.
                currentSubstring += currChar;
            }
            else if(currChar == CHAR::TIMES && stringSpot > 0 && _functionString[stringSpot-1] == CHAR::TIMES) {
                //Do nothing, this is a ** which is a ^.
                currentSubstring += currChar;
            }
            //Store the current substring and the new sign if it exists
            else if(currChar == CHAR::TIMES || currChar == CHAR::DIVIDE) {
                _subparts.push_back(currentSubstring);
                _isMultiply.push_back(currentIsMultiply);
                currentSubstring = "";
                currentIsMultiply = currChar == CHAR::TIMES;
            }
            else {
                currentSubstring += currChar;
            }
        }
        
        //Add Anything Remaining
        _subparts.push_back(currentSubstring);
        _isMultiply.push_back(currentIsMultiply);
    }

    void splitPower(const std::string& _functionString, std::vector<std::string>& _subparts) {
        std::string currentSubstring = "";
        size_t parenthesisCount = 0;
        for (std::string::size_type stringSpot = 0; stringSpot < _functionString.size(); stringSpot++) {
            char currChar = _functionString[stringSpot];
            //Update the parenthesisCount
            if(currChar == CHAR::RIGHT_PAREN) {
                if(parenthesisCount == 0) {
                    printAndThrowRuntimeError("Mismatched Parenthesis in Function!");
                }
                parenthesisCount--;
                currentSubstring += currChar;
            }
            else if(currChar == CHAR::LEFT_PAREN) {
                parenthesisCount++;
                currentSubstring += currChar;
            }
            else if(parenthesisCount > 0) {
                //Do nothing, we are inside parenthesis
                currentSubstring += currChar;

            }
            //Check for a power split ^
            else if(currChar == CHAR::POWER) {
                _subparts.push_back(currentSubstring);
                currentSubstring = "";
            }
            //Check for a power split **
            else if(currChar == CHAR::TIMES && stringSpot + 1 < _functionString.length() && _functionString[stringSpot+1] == CHAR::TIMES) {
                _subparts.push_back(currentSubstring);
                currentSubstring = "";
                stringSpot++; //Skip past the second CHAR::TIMES
            }
            else {
                currentSubstring += currChar;
            }
        }
        
        //Add Anything Remaining
        _subparts.push_back(currentSubstring);
    }

    void parseSum(const std::vector<std::string>& _subparts) {
        m_functionType = FunctionType::SUM;
        m_value = 0;
        //Add all the subfunctions
        for(size_t i = 0; i < _subparts.size(); i++) {
            addSubfunction(_subparts[i]);
        }
    }

    void parseProduct(const std::vector<std::string>& _subparts) {
        m_functionType = FunctionType::PRODUCT;
        m_value = 1;
        //Add all the subfunctions
        for(size_t i = 0; i < _subparts.size(); i++) {
            addSubfunction(_subparts[i]);
        }
    }

    void parsePower(const std::vector<std::string>& _subparts) {
        //Power must have exactly 2 subparts!
        if(_subparts.size() != 2) {
            printAndThrowRuntimeError("Error Parsing Functions! More than 1 power sign used in string! Note: 2^2^2 is not allowed, split up with parenthesis!");
        }
        
        m_functionType = FunctionType::POWER;
        m_value = 1;
        //Add all the subfunctions
        for(size_t i = 0; i < _subparts.size(); i++) {
            addSubfunction(_subparts[i]);
        }
    }

    bool parseNumber(std::string& _functionString) {
        //Checks if this is a valid number.
        try{
            //TODO: std::stod doesn't always throw and error if things are wrong.
            //      Example: 2(x-1) parses as 2.
            //      Write my own parser that will parse from numbers to doubles!
            m_value = std::stod(_functionString);
            return true;
        }
        catch(...) {
            return false;
        }
    }
    
    void parseComplexType(std::string& _functionString) {
        //These all start being multiplied by 1. They have the ability to have an m_value multiplied by it so
        //in simplifyExpressions we could pull constants in the multiply out and put it in m_value.
        //TODO: Is this a good idea or should I just drop m_value from complex things?
        m_value = 1;
        
        //Parse constants
        if(parseNumber(_functionString)) { //Numbers
            m_functionType = FunctionType::CONSTANT;
            return;
        }
        if(_functionString == "e" || _functionString == "E") { //e
            m_functionType = FunctionType::CONSTANT;
            m_value *= M_E;
            return;
        }
        if(_functionString == "pi" || _functionString == "PI" || _functionString == "Pi" || _functionString == "pI") { //pi
            m_functionType = FunctionType::CONSTANT;
            m_value *= M_PI;
            return;
        }
        
        size_t stringLenth = _functionString.length();
        
        //Parse 5 char subfunctions: log10
        if(stringLenth > 5) {
            std::string lowercase5 = toLowerSubstring(_functionString, 0, 5);
            if(lowercase5 == "log10") {
                m_functionType = FunctionType::LOG10;
                parseComplexFunctionParenthesis(_functionString.substr(5));
                return;
            }
        }
        
        //Parse 4 char subfunctions: hyperbolic trig, sqrt, log2 and prod
        if(stringLenth > 4) {
            std::string lowercase4 = toLowerSubstring(_functionString, 0, 4);
            if(lowercase4 == "sinh") {
                m_functionType = FunctionType::SINH;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "cosh") {
                m_functionType = FunctionType::COSH;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "tanh") {
                m_functionType = FunctionType::TANH;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "sqrt") {
                m_functionType = FunctionType::SQRT;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "log2") {
                m_functionType = FunctionType::LOG2;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "prod") {
                m_value = 1;
                m_functionType = FunctionType::PRODUCT;
                parseComplexFunctionExpansionParenthesis(_functionString.substr(4));
                return;
            }
        }
        
        //Parse 3 char subfunctions: Normal Trig, exp, log and sum
        if(stringLenth > 3) {
            std::string lowercase3 = toLowerSubstring(_functionString, 0, 3);
            if(lowercase3 == "sin") {
                m_functionType = FunctionType::SIN;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "cos") {
                m_functionType = FunctionType::COS;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "tan") {
                m_functionType = FunctionType::TAN;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "exp") {
                m_functionType = FunctionType::EXP;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "log") {
                m_functionType = FunctionType::LOG;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "sum") {
                m_value = 0;
                m_functionType = FunctionType::SUM;
                parseComplexFunctionExpansionParenthesis(_functionString.substr(3));
                return;
            }
        }
        
        //Parse VARIABLE
        for(size_t i = 0; i < m_variableNames.size(); i++) {
            if(m_variableNames[i] == _functionString) {
                m_functionType = FunctionType::VARIABLE;
                m_varIndex = i;
                return;
            }
        }
        
        //Parse CHEBYSHEV
        if(_functionString.length() > 1 && std::tolower(_functionString[0]) == 't' && isNumericDigit(_functionString[1])) {
            parseChebyshevFunction(_functionString);
            m_functionType = FunctionType::CHEBYSHEV;
            return;
        }
                
        //It's not a function we recogize, throw an error!
        printAndThrowRuntimeError("Failed to Parse Function! Unknown Type: " + _functionString);
    }
    
    void parseComplexFunctionParenthesis(const std::string& _functionString) {
        //Following most complex functions is (<subfunction>). Strip the () and add the subfunction.
        if(_functionString[0] != CHAR::LEFT_PAREN || _functionString[_functionString.length() - 1] != CHAR::RIGHT_PAREN) {
            printAndThrowRuntimeError("Failed to Parse Function!");
        }
        addSubfunction(_functionString.substr(1, _functionString.length() - 2));
    }
    
    void parseComplexFunctionExpansionParenthesis(const std::string& _functionString) {
        //This is for sum and prod expansions. The are of the form <sum/prod>(<function>, <var>, <start>, <end>)
        
        //Make sure the parenthesis are there.
        if(_functionString[0] != CHAR::LEFT_PAREN || _functionString[_functionString.length() - 1] != CHAR::RIGHT_PAREN) {
            printAndThrowRuntimeError("Failed to Parse Function!");
        }
        //Split up the subparts on the comma
        std::vector<std::string> parts = split(_functionString.substr(1, _functionString.length() - 2), ",");
        if(parts.size() != 4) {
            printAndThrowRuntimeError("Failed to Parse Function! Incorrect format!");
        }
        //Make sure the <var> isn't also the name of some function, constant, or variable, as the causes ambiguity.
        checkNameNotClaimed(parts[1], "Sum/Prod Var");
        
        //Get the subtrings we will use.
        std::vector<std::string> subFunctionStrings = split(parts[0], parts[1]);
        //If subFunctionStrings[i][-1] and subFunctionStrings[i+1][0] are not both ()+-*/^ types, then combine them.
        //This avoids replacing the i in sin for example, or other variable names.
        size_t indexToCheck = 1;
        while(indexToCheck < subFunctionStrings.size()) { //Note that subFunctionStrings could be size 0 at the end of this.
            //Make sure the subFunctionStrings have size
            if(subFunctionStrings[indexToCheck-1].length() == 0) {
                subFunctionStrings.erase(subFunctionStrings.begin() + (indexToCheck-1));
                continue;
            }
            else if(subFunctionStrings[indexToCheck].length() == 0) {
                subFunctionStrings.erase(subFunctionStrings.begin() + indexToCheck);
                continue;
            }
            
            //Check if the chars are valid delimiters
            char lastOfFirst = subFunctionStrings[indexToCheck-1][subFunctionStrings[indexToCheck-1].length() - 1];
            char firstOfLast = subFunctionStrings[indexToCheck][0];
            const bool firstGood = lastOfFirst == CHAR::LEFT_PAREN || lastOfFirst == CHAR::RIGHT_PAREN || lastOfFirst == CHAR::PLUS || lastOfFirst == CHAR::MINUS || lastOfFirst == CHAR::DIVIDE || lastOfFirst == CHAR::TIMES || lastOfFirst == CHAR::POWER;
            const bool lastGood = firstOfLast == CHAR::LEFT_PAREN || firstOfLast == CHAR::RIGHT_PAREN || firstOfLast == CHAR::PLUS || firstOfLast == CHAR::MINUS || firstOfLast == CHAR::DIVIDE || firstOfLast == CHAR::TIMES || firstOfLast == CHAR::POWER;
            //Recombine if needed, otherwise move on.
            if(!firstGood || !lastGood) {
                subFunctionStrings[indexToCheck-1] += parts[1] + subFunctionStrings[indexToCheck];
                subFunctionStrings.erase(subFunctionStrings.begin() + indexToCheck);
            }
            else {
                indexToCheck++;
            }
        }
        
        //Get the bounds for the sum/prod
        int64_t start, end;
        try {
            start = std::stoll(parts[2]);
            end = std::stoll(parts[3]);
        }
        catch(...) {
            printAndThrowRuntimeError("Failed to Parse Function! Illegal Bound Value in sum/prod, not an integer");
        }
        
        //Add the subfunctions with the number replaced
        m_operatorSigns.clear();
        for(int64_t replaceNum = start; replaceNum <= end; replaceNum++) {
            std::string resultFunc = subFunctionStrings[0];
            std::string replaceString = "(" + std::to_string(replaceNum) + ")";
            for(size_t i = 1; i < subFunctionStrings.size(); i++) {
                resultFunc += replaceString + subFunctionStrings[i];
            }
            addSubfunction(resultFunc);
            m_operatorSigns.push_back(true); //This is always + or *.
        }
    }

    void parseChebyshevFunction(const std::string& _functionString) { //TODO: This could be updated and needs better comments.
        if(_functionString[_functionString.length() - 1] != CHAR::RIGHT_PAREN) {
            printAndThrowRuntimeError("Failed to Parse Function!");
        }
        std::string substring1;
        std::string substring2;
        bool part2 = false;
        size_t parenthesisCount = 0;
        //The two subfunctions are seperated by a comma
        for(size_t i = 1; i + 1 < _functionString.length(); i++) {
            char c = _functionString[i];
            //Update the parenthesisCount
            if(c == CHAR::RIGHT_PAREN) {
                if(parenthesisCount == 0) {
                    printAndThrowRuntimeError("Mismatched Parenthesis in Function!");
                }
                parenthesisCount--;
            }
            else if(c == CHAR::LEFT_PAREN && part2) {
                parenthesisCount++;
            }
            //Check for a comma and add to the right substring
            if(c == CHAR::LEFT_PAREN && !part2) {
                part2 = true;
            }
            else if (part2) {
                substring2 += c;
            }
            else {
                substring1 += c;
            }
        }
        
        m_varIndex = std::stoull(substring1);
        addSubfunction(substring2);
    }

    bool parsePower(const std::string& _functionString, const std::string& _coeffString) {
        std::string substring1;
        std::string substring2;
        bool part2 = false;
        size_t parenthesisCount = 0;
        //The two subfunctions are seperated by a comma
        for(size_t i = 0; i < _functionString.length(); i++) {
            char c = _functionString[i];
            //Update the parenthesisCount
            if(c == CHAR::RIGHT_PAREN) {
                if(parenthesisCount == 0) {
                    printAndThrowRuntimeError("Mismatched Parenthesis in Function!");
                }
                parenthesisCount--;
            }
            else if(c == CHAR::LEFT_PAREN) {
                parenthesisCount++;
            }
            //Check for a power symbol and split the substrings
            if(parenthesisCount == 0 && c == CHAR::POWER) {
                part2 = true;
            }
            else if(parenthesisCount == 0 && c == CHAR::TIMES && i+1 < _functionString.length() && _functionString[i+1] == CHAR::TIMES) {
                //** is parsed as ^
                part2 = true;
                i++;
            }
            else if (part2) {
                substring2 += c;
            }
            else {
                substring1 += c;
            }
        }
        if(part2) {
            if(substring1 == "") {//The coeff might have been the whole rhs, so put it back.
                //If the _coeffString has a - in it, pull that out and make the m_value -1. Because then negative should not be in the power.
                if(_coeffString[0] == CHAR::MINUS) {
                    substring1 = _coeffString.substr(1);
                    m_value = -1;
                }
                else {
                    substring1 = _coeffString;
                    m_value = 1;
                }
            }
            addSubfunction(substring1);
            addSubfunction(substring2);
            return true;
        }
        return false;
    }
    
    void addSubfunction(const std::string& subfunctionString) {
        if(unlikely(s_allFunctions.size() == 0)) {
            s_allFunctions.resize(1);
            s_allFunctionNames.resize(1);
        }

        //Get rid of a + at the start of the string
        std::string newString = subfunctionString;
        if(newString.length() > 0 && newString[0] == CHAR::PLUS) {
            newString = newString.substr(1);
        }
        
        //Remove parenthesis
        removeExtraParenthesis(newString);
        
        //Check if this is already a function
        FunctionMap::const_iterator found = s_allFunctions[0].find(newString);
        if(found != s_allFunctions[0].end()) {
            m_subfunctions.push_back(found->second);
            return;
        }
        //Check if this is already a function name.
        found = s_allFunctionNames[0].find(newString);
        if(found != s_allFunctionNames[0].end()) {
            m_subfunctions.push_back(found->second);
            return;
        }
        
        //Make a new function
        m_subfunctions.push_back(addFunction("", newString, m_variableNames));
    }
    
//Specialized Function Evals
    template<typename ReturnType>
    ReturnType sumEval(const std::vector<double>& _inputPoints) {
        assert(m_subfunctions.size() == m_operatorSigns.size());
        ReturnType result = m_value;
        for(size_t i = 0; i < m_subfunctions.size(); i++) {
            SIGNCHECKSUM(m_operatorSigns[i], result, m_subfunctions[i]->evaluate<ReturnType>(_inputPoints));
        }
        return result;
        
        //TODO: Add stable eval option where a sum eval gets all the things to sum and then does it from smallest to largest.
        //Maybe even make this the default?
        //Also, would it be possible to track the potential error inside the evaluation?
        //Like when adding two numbers we know the error is machine epsilon * the larger number.
        //Have a STABLE_SUM type that does stable sum evals
        //Possibly instead of summing doubles sum type that is templated. This can be doubles, or can be a struct that
        //is a double and then an error, so the end result of the eval returns the answer and the maximum error.
    }
    
    template<typename ReturnType>
    ReturnType productEval(const std::vector<double>& _inputPoints) {
        assert(m_subfunctions.size() == m_operatorSigns.size());
        ReturnType result = m_value;
        for(size_t i = 0; i < m_operatorSigns.size(); i++) {
            SIGNCHECKPRODUCT(m_operatorSigns[i], result, m_subfunctions[i]->evaluate<ReturnType>(_inputPoints));
        }
        return result;
    }
    
public:
    void getFunctionLevels(std::vector<std::vector<SharedFunctionPtr> >& allFunctionLevels) {
        //The function level is 1 + max(function level of subfunctions).
        //If a function has no subfunctions, the function level is 0.
        //This way when we evaluate we can evaluate all the functions starting a level 0 at going up.
        
        m_functionLevel = 0;
        for(size_t subFuncNum = 0; subFuncNum < m_subfunctions.size(); subFuncNum++) {
            m_subfunctions[subFuncNum]->getFunctionLevels(allFunctionLevels);
            size_t currLevel = m_subfunctions[subFuncNum]->getFunctionLevel();
            m_functionLevel = std::max(m_functionLevel, currLevel + 1);
            
            if(allFunctionLevels.size() <= currLevel) {
                allFunctionLevels.resize(currLevel + 1);
            }
            
            //Add to the vector if it's not already there
            bool alreadyExists = false;
            for(size_t i = 0; i < allFunctionLevels[currLevel].size(); i++) {
                if(allFunctionLevels[currLevel][i] == m_subfunctions[subFuncNum]) {
                    alreadyExists = true;
                    break;
                }
            }
            if(!alreadyExists) {
                allFunctionLevels[currLevel].push_back(m_subfunctions[subFuncNum]);
            }
            
        }
        
        //Figure out what dimensions the function exists in.
        defineFunctionDimensions();
    }
        
//Overloaded Operators
public:
    friend inline bool operator == (const Function& lhs, const Function& rhs){ return lhs.getFunctionString() == rhs.getFunctionString(); }
    friend inline bool operator != (const Function& lhs, const Function& rhs){ return !(lhs == rhs); }
    
    friend inline bool operator< (const Function& lhs, const Function& rhs){ return lhs.getFunctionString() < rhs.getFunctionString(); }
    friend inline bool operator> (const Function& lhs, const Function& rhs){ return rhs < lhs; }
    friend inline bool operator<=(const Function& lhs, const Function& rhs){ return !(lhs > rhs); }
    friend inline bool operator>=(const Function& lhs, const Function& rhs){ return !(lhs < rhs); }
    
//Getters
public:
    const std::vector<bool>& getOperatorSigns() const {
        return m_operatorSigns;
    }

    FunctionType getFunctionType() const {
        return m_functionType;
    }
    
    double getValue() const {
        return m_value;
    }
    
    size_t getVarIndex() const {
        return m_varIndex;
    }
    
    const std::string& getFunctionString() const {
        return m_functionString;
    }

    const std::string& getFunctionName() const {
        return m_functionName;
    }

    const std::vector<std::string>& getVariableNames() const {
        return m_variableNames;
    }

    const Monomial& getMonomial() const {
        return m_monomial;
    }

    const Polynomial& getPolynomial() const {
        return m_polynomial;
    }
    
    bool isTopFunction() const {
        return m_isTopFunction;
    }
    
    size_t getFunctionLevel() const {
        return m_functionLevel;
    }
    
    const std::vector<double>& getPartialEvals() const {
        return m_partialEvaluations;
    }

    const std::vector<bool>& getHasDimension() const {
        return m_hasDimension;
    }

    size_t getRank() const {
        return m_rank;
    }

    size_t getNumUsedDimensions() const {
        return m_numUsedDimensions;
    }

    const std::vector<EvaluateGridInfo>& getEvaluateGridInfo() const {
        return m_evaluateGridInfo;
    }

    EvaluateGridType getEvaluateGridType() const {
        return m_evaluateGridType;
    }

    const std::vector<std::vector<SharedFunctionPtr> >& getAllFunctionLevels() const {
        return m_allFunctionLevels;
    }

    const std::vector<SharedFunctionPtr>& getSubfunctions() const {
        return m_subfunctions;
    }
        
    static SharedFunctionPtr getThreadFunctionByName(size_t _threadNum, const std::string& _functionName) {
        if(_threadNum >= s_allFunctionNames.size()) {
            printAndThrowRuntimeError("Functions not copied for thread number " + std::to_string(_threadNum));
        }
        if(s_allFunctionNames[_threadNum].find(_functionName) == s_allFunctionNames[_threadNum].end()) {
            printAndThrowRuntimeError("No definition found for function " + _functionName + " on thread " + std::to_string(_threadNum));
        }
        return s_allFunctionNames[_threadNum][_functionName];
    }

    static void addThreadFunctions(size_t _requiredThreads) {
        while(s_allFunctions.size() <_requiredThreads) {
            //Add 1 to the vectors
            s_allFunctions.resize(s_allFunctions.size() + 1);
            s_allFunctionNames.resize(s_allFunctionNames.size() + 1);
            const size_t size = s_allFunctions.size()-1;
            
            //Copy all s_allFunctions.
            for(FunctionMap::const_iterator it = s_allFunctions[0].begin(); it != s_allFunctions[0].end(); it++) {
                //Copy this function
                SharedFunctionPtr currFunction = std::make_shared<Function>(*it->second.get());
                //Add it to s_allFunctions
                s_allFunctions[size].insert({it->first,currFunction});
            }
            
            //Finalize Copy on Everything
            for(FunctionMap::const_iterator it = s_allFunctions[size].begin(); it != s_allFunctions[size].end(); it++) {
                //Find the function to copy by string, and copy it
                std::string funcStringToFind = it->first;
                FunctionMap::const_iterator found = s_allFunctions[0].find(funcStringToFind);
                std::shared_ptr<Function> functionToCopy = found->second;
                it->second->finalizeCopy(*functionToCopy.get());
            }
            
            //Copy s_allFunctionNames
            for(FunctionMap::const_iterator it = s_allFunctions[size].begin(); it != s_allFunctions[size].end(); it++) {
                const std::string& functionName = it->second->getFunctionName();
                if(functionName != "") { //The should be every function that has a name.
                    s_allFunctionNames[s_allFunctionNames.size() - 1][functionName] = it->second;
                }
            }
        }
    }
    
    //Copy Constructor that doesn't copy the pointers.
    Function(const Function& other) :
    m_isTopFunction(other.isTopFunction()),
    m_functionLevel(other.getFunctionLevel()),
    m_operatorSigns(other.getOperatorSigns()),
    m_functionType(other.getFunctionType()),
    m_value(other.getValue()),
    m_varIndex(other.getVarIndex()),
    m_functionString(other.getFunctionString()),
    m_functionName(other.getFunctionName()),
    m_variableNames(other.getVariableNames()),
    m_monomial(other.getMonomial()),
    m_polynomial(other.getPolynomial()),
    m_rank(other.getRank()),
    m_hasDimension(other.getHasDimension()),
    m_numUsedDimensions(other.getNumUsedDimensions()),
    m_partialEvaluations(other.getPartialEvals()),
    m_evaluateGridInfo(other.getEvaluateGridInfo()),
    m_evaluateGridType(other.getEvaluateGridType())
    {}
    
    //Helper to the copy constructor, updates the pointers
    void finalizeCopy(const Function& other) {
        const size_t threadNum = s_allFunctions.size() - 1; //Copy from this spot in the s_allFunctions
        
        //Copy m_allFunctionLevels
        const std::vector<std::vector<SharedFunctionPtr> >& functionLevelsToCopy = other.getAllFunctionLevels();
        m_allFunctionLevels.resize(functionLevelsToCopy.size());
        for(size_t level = 0; level < functionLevelsToCopy.size(); level++) {
            m_allFunctionLevels[level].resize(functionLevelsToCopy[level].size());
            for(size_t innerLevel = 0; innerLevel < functionLevelsToCopy[level].size(); innerLevel++) {
                m_allFunctionLevels[level][innerLevel] = s_allFunctions[threadNum][functionLevelsToCopy[level][innerLevel]->getFunctionString()];
            }
        }
        
        //Copy m_subfunctions
        const std::vector<SharedFunctionPtr>& subfunctionsToCopy = other.getSubfunctions();
        m_subfunctions.resize(subfunctionsToCopy.size());
        for(size_t level = 0; level < subfunctionsToCopy.size(); level++) {
            m_subfunctions[level] = s_allFunctions[threadNum][subfunctionsToCopy[level]->getFunctionString()];
        }
    }
    
private:
    static void checkNameNotClaimed(const std::string& _name, const std::string& _type) {
        if(s_claimedVariableNames.find(_name) != s_claimedVariableNames.end() && _type != "Variable Name") {
            printAndThrowRuntimeError("Illegal " + _type + " Name: " + _name + ". Already a Variable Name.");
        }
        else if(s_allFunctionNames[0].find(_name) != s_allFunctionNames[0].end()) {
            printAndThrowRuntimeError("Illegal " + _type + " Name: " + _name + ". Already a Function Name.");
        }
        else if(s_claimedFunctionNames.find(_name) != s_claimedFunctionNames.end()) {
            printAndThrowRuntimeError("Illegal " + _type + " Name: " + _name + ". Already a Set Function Name.");
        }
        else if(s_claimedConstantNames.find(_name) != s_claimedConstantNames.end()) {
            printAndThrowRuntimeError("Illegal " + _type + " Name: " + _name + ". Already a Set Constant Name.");
        }
    }
    
public:
    //A function should only be created through this. It adds them to the maps correctly.
    static SharedFunctionPtr addFunction(const std::string& _functionName, const std::string& _functionString, const std::vector<std::string>& _variableNames) {
        if(unlikely(s_allFunctions.size() == 0)) {
            s_allFunctions.resize(1);
            s_allFunctionNames.resize(1);
        }
        
        //Claim the variable names
        for(size_t i = 0; i < _variableNames.size(); i++) {
            checkNameNotClaimed( _variableNames[i], "Variable Name");
            s_claimedVariableNames.insert(_variableNames[i]);
        }
        
        //Make sure this is a valid name.
        checkNameNotClaimed(_functionName, "Function");

        //Check if we are double adding it
        bool namedFuncAlreadyExists = false;
        if(_functionName != "") {
            if(s_allFunctions[0].find(_functionString) != s_allFunctions[0].end()) {
                namedFuncAlreadyExists = true;
            }
        }
        else if(s_allFunctions[0].find(_functionString) != s_allFunctions[0].end()) {
            printAndThrowRuntimeError("Trying to add function Twice! Name is " + _functionName + " String is " + _functionString);
        }
        
        if(namedFuncAlreadyExists) {
            SharedFunctionPtr function = s_allFunctions[0][_functionString];
            s_allFunctionNames[0][_functionName] = function;
            return function;
        }
        else {
            SharedFunctionPtr function = std::make_shared<Function>(_functionName, _functionString, _variableNames);

            if(_functionName != "") {
                s_allFunctionNames[0][_functionName] = function;
            }
            s_allFunctions[0][_functionString] = function;
            
            return function;
        }
    }
        
    static void clearSavedFunctions() {
        s_allFunctions.clear();
        s_allFunctionNames.clear();
    }

private:
    //The Top Function and m_allFunctions Map.
    static std::vector<FunctionMap>             s_allFunctions; //Map from function string to function
    static std::vector<FunctionMap>             s_allFunctionNames; //Map from function name to function
    
    bool                                        m_isTopFunction;
    size_t                                      m_functionLevel;
    std::vector<std::vector<SharedFunctionPtr> > m_allFunctionLevels;

    //The subfunctions and coresponding signs
    std::vector<SharedFunctionPtr>              m_subfunctions;
    std::vector<bool>                           m_operatorSigns; //True is + or *. False is - or /.
    
    //The type of function
    FunctionType                                m_functionType;
    
    //The Signed Coefficient of anything that isn't a sum. Sums all have coefficient of 1.
    double                                      m_value;
    
    //For VARIABLE: The variable index
    size_t                                      m_varIndex;
            
    //Variable Names and Subfunctions Names
    std::string                                 m_functionString;
    std::string                                 m_functionName;
    std::vector<std::string>                    m_variableNames;
    
    //For Polynomials
    Monomial                                    m_monomial;
    Polynomial                                  m_polynomial;
    
    //For evaluate grid
    size_t                                      m_rank;
    std::vector<bool>                           m_hasDimension;
    size_t                                      m_numUsedDimensions;
    std::vector<double>                         m_partialEvaluations;
    std::vector<EvaluateGridInfo>               m_evaluateGridInfo;
    EvaluateGridType                            m_evaluateGridType;
        
    static                  std::unordered_set<std::string> s_claimedConstantNames;
    static                  std::unordered_set<std::string> s_claimedFunctionNames;
    static                  std::unordered_set<std::string> s_claimedVariableNames;
};

std::vector<Function::FunctionMap>  Function::s_allFunctions;
std::vector<Function::FunctionMap>  Function::s_allFunctionNames;

std::unordered_set<std::string> Function::s_claimedConstantNames({"e","pi",});
std::unordered_set<std::string> Function::s_claimedFunctionNames({"sin","cos","tan","sinh","cosh","tanh","sqrt","exp","log","log2","log10"});
std::unordered_set<std::string> Function::s_claimedVariableNames;


#endif /* Function_h */
