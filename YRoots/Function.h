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
#include <unordered_map>
#include <unordered_set>
#include "Timer.h"

//TODO: Maybe I should set this up so that coeffs are stored seperately on sums and products? So 4cos(x)+5cos(x) use the same cos(x)?

enum FunctionType {
    SIN, //Syntax: sin(x)
    COS, //Syntax: cos(x)
    TAN, //Syntax: tan(x)
    COSH, //Syntax: cosh(x)
    SINH, //Syntax: sinh(x)
    TANH, //Syntax: tanh(x)
    LOG, //Syntax: log(base, x)
    LN, //Syntax: ln(x)
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
    std::vector<std::vector<size_t>> childEvalIndexes;
    size_t childEvalSize;
};

#define SIGNCHECKPRODUCT(isPositive, sum, number) (isPositive ? sum *= number : sum /= number)

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
    m_dimension(_variableNames.size())
    {
        //Parse the function
        functionParse(m_functionString);
        
        m_timer.registerTimer(m_timerEvaluateGridIndex, "Evaluate Grid");
    }
    
    void defineFunctionDimensions() {
        //Populate m_hasDimension
        m_hasDimension.resize(m_dimension);
        switch(m_functionType) {
            case FunctionType::CONSTANT:
                break;
            case FunctionType::VARIABLE:
                m_hasDimension[m_varIndex] = true;
                break;
            case FunctionType::POWER_BASIS_MONOMIAL: //TODO: Write this for monomials and polynomials
                break;
            case FunctionType::CHEBYSHEV_BASIS_MONOMIAL:
                break;
            case FunctionType::POWER_BASIS_POLYNOMIAL:
                break;
            case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL:
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
            case FunctionType::POWER_BASIS_MONOMIAL: //TODO: Write the specialized base calls for monomials and polynomials
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
            case FunctionType::LN:
            case FunctionType::SQRT:
            case FunctionType::EXP:
            case FunctionType::CHEBYSHEV:
                m_evaluateGridType = EvaluateGridType::SIMPLE;
                break;
            case FunctionType::POWER:
            case FunctionType::LOG:
            case FunctionType::SUM:
            case FunctionType::PRODUCT:
                m_evaluateGridType = EvaluateGridType::COMBINE;
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid main! Fix Switch Statement!");
                break;
        }
        
        //The top function can only be simple if it has every dimension
        if(m_isTopFunction && m_evaluateGridType == EvaluateGridType::SIMPLE) {
            for(size_t i = 0; i < m_hasDimension.size(); i++) {
                if(!m_hasDimension[i]) {
                    m_evaluateGridType = EvaluateGridType::COMBINE;
                    break;
                }
            }
        }
    }
    
    void prepEvaluatGrid(size_t _gridSize, EvaluateGridInfo& _infoToPopulate) {
        _infoToPopulate.precomputed = true;
        
        //Both cases we need to know the childEvalSize
        _infoToPopulate.childEvalSize = power(_gridSize, m_numUsedDimensions);

        //Make sure the m_partialEvaluations is big enough
        if(m_partialEvaluations.size() <= _infoToPopulate.childEvalSize) {
            m_partialEvaluations.resize(_infoToPopulate.childEvalSize);
        }

        if(m_evaluateGridType == EvaluateGridType::COMBINE) {
            //We need to create the childEvalIndexes for the non simple cases
            _infoToPopulate.childEvalIndexes.resize(_infoToPopulate.childEvalSize);
            
            //Start on the first dimensions
            std::vector<size_t> tempDimCount(m_dimension);
            size_t dimToInc = 0;
            if(!m_isTopFunction){
                while(!m_hasDimension[dimToInc] && dimToInc < m_dimension) {
                    dimToInc++;
                }
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
                    //Go the the nest spot
                    if(!m_isTopFunction){
                        while(!m_hasDimension[dimToInc] && dimToInc < m_dimension) {
                            dimToInc++;
                        }
                    }
                    //Check if it overflows as well
                    if(dimToInc < m_dimension) {
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
    }
    
    void evaluateGridBase(const std::vector<std::vector<double>>& _grid) {
        switch(m_functionType) {
            case FunctionType::CONSTANT:
                m_partialEvaluations[0] =  m_value;
                break;
            case FunctionType::VARIABLE:
                for(size_t i = 0; i < _grid[0].size(); i++) {
                    m_partialEvaluations[i] =  m_value * _grid[m_varIndex][i];
                }
                break;
            case FunctionType::POWER_BASIS_MONOMIAL: //TODO: Write the specialized base calls for monomials and polynomials
                break;
            case FunctionType::CHEBYSHEV_BASIS_MONOMIAL:
                break;
            case FunctionType::POWER_BASIS_POLYNOMIAL:
                break;
            case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL:
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid base! Fix Switch Statement!");
                break;
        }
    }

    void evaluateGridBaseTop(const std::vector<std::vector<double>>& _grid, std::vector<double>& _results) {
        size_t gridSize = _grid[0].size();
        size_t numEvals = power(gridSize, _grid.size());;
        
        switch(m_functionType) {
            case FunctionType::CONSTANT:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value;
                }
                break;
            case FunctionType::VARIABLE: {
                //TODO: this could probably be more efficient.
                size_t varGridSize = power(gridSize, m_varIndex);
                for(size_t i = 0; i < numEvals; i++) {
                    size_t index = (i/varGridSize)%gridSize;
                    _results[i] =  _grid[m_varIndex][index];
                }
                break;
            }
            case FunctionType::POWER_BASIS_MONOMIAL: //TODO: Write the specialized base calls for monomials and polynomials
                break;
            case FunctionType::CHEBYSHEV_BASIS_MONOMIAL:
                break;
            case FunctionType::POWER_BASIS_POLYNOMIAL:
                break;
            case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL:
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid base! Fix Switch Statement!");
                break;
        }
    }

    void evaluateGridSimple(size_t _gridSize, std::vector<double>& _results) {
        //Get the childs evaluations and resize the result vector if needed
        const std::vector<double>& childEvals = m_subfunctions[0]->getPartialEvals();
        size_t numEvals = m_evaluteGridInfo[_gridSize].childEvalSize;

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
            case FunctionType::LN:
                for(size_t i = 0; i < numEvals; i++) {
                    _results[i] =  m_value * log(childEvals[i]);
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
        std::vector<std::vector<size_t>>& currSpots = m_evaluteGridInfo[_gridSize].childEvalIndexes;
        switch(m_functionType) {
            case FunctionType::POWER: {
                const std::vector<double>& childEvals1 = m_subfunctions[0]->getPartialEvals();
                const std::vector<double>& childEvals2 = m_subfunctions[1]->getPartialEvals();
                for(size_t i = 0; i < currSpots.size(); i++) {
                    _results[i] =  m_value * pow(childEvals1[currSpots[i][0]], childEvals2[currSpots[i][1]]);
                }
                break;
            }
            case FunctionType::LOG:{
                const std::vector<double>& childEvals1 = m_subfunctions[0]->getPartialEvals();
                const std::vector<double>& childEvals2 = m_subfunctions[1]->getPartialEvals();
                for(size_t i = 0; i < currSpots.size(); i++) {
                    //TOOD: This could be more efficent by just taking the logs for each dimensions, and then doing the division
                    //In the cross product, but that's not a priority as the function is only used for things of variable base,
                    //which I expect to be rare. Maybe I should just turn things of variable base into a product type and then
                    //it isn't an issue.
                    _results[i] =  m_value * log(childEvals1[currSpots[i][0]]) / log(childEvals2[currSpots[i][1]]);
                }
                break;
            }
            case FunctionType::SUM:{
                double result;
                for(size_t i = 0; i < currSpots.size(); i++) {
                    result = m_value;
                    for(size_t j = 0; j < currSpots[i].size(); j++) {
                        result += m_subfunctions[j]->getPartialEvals()[currSpots[i][j]];
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
                        SIGNCHECKPRODUCT(m_isMultiply[j], result, m_subfunctions[j]->getPartialEvals()[currSpots[i][j]]);
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
    
    void evaluateGridMain(const std::vector<std::vector<double>>& _grid) {
        //Prep the evaluate grid info if it hasn't already been done
        size_t gridSize = _grid[0].size();
        if(unlikely(m_evaluteGridInfo.size() <= gridSize)) {
            m_evaluteGridInfo.resize(gridSize + 1);
        }
        if(unlikely(!m_evaluteGridInfo[gridSize].precomputed)) {
            prepEvaluatGrid(gridSize, m_evaluteGridInfo[gridSize]);
        }

        switch(m_evaluateGridType) {
            case EvaluateGridType::BASE:
                evaluateGridBase(_grid);
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

    void evaluateGridMain(const std::vector<std::vector<double>>& _grid, std::vector<double>& _results) {
        //Prep the evaluate grid info if it hasn't already been done
        size_t gridSize = _grid[0].size();
        if(unlikely(m_evaluteGridInfo.size() <= gridSize)) {
            m_evaluteGridInfo.resize(gridSize + 1);
        }
        if(unlikely(!m_evaluteGridInfo[gridSize].precomputed)) {
            prepEvaluatGrid(gridSize, m_evaluteGridInfo[gridSize]);
        }

        switch(m_evaluateGridType) {
            case EvaluateGridType::BASE:
                evaluateGridBaseTop(_grid, _results);
                break;
            case EvaluateGridType::SIMPLE:
                evaluateGridSimple(gridSize, _results);
                break;
            case EvaluateGridType::COMBINE:
                evaluateGridCombine(gridSize, _results);
                break;
            default:
                printAndThrowRuntimeError("Unknown Function Type Encountered in evaluate grid main! Fix Switch Statement!");
                break;
        }
    }
    
    void evaluateGrid(const std::vector<std::vector<double>>& _grid, std::vector<double>& _results) {
        //Evaluates the grid so that the evaluation of (grid[0][i], grid[1][j]) is in _results[j+grid.size() + i]
        m_timer.startTimer(m_timerEvaluateGridIndex);
        
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
        m_timer.stopTimer(m_timerEvaluateGridIndex);
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
        //assert(_results[evalSpot++] == evaluate(inputPoints));
        _results[evalSpot++] = evaluate(inputPoints);
        while (spotToInc < dimension) {
            bool firstPass = true;
            while(++inputSpot[spotToInc] < numPoints) {
                inputPoints[spotToInc] = _grid[spotToInc][inputSpot[spotToInc]];
                //assert(_results[evalSpot++] == evaluate(inputPoints));
                _results[evalSpot++] = evaluate(inputPoints);
                if(firstPass && spotToInc != 0) {
                    spotToInc = 0;
                }
                firstPass = false;
            }
            inputSpot[spotToInc] = 0;
            inputPoints[spotToInc] = _grid[spotToInc][0];
            spotToInc++;
        }
    }
    
    double evaluate(const std::vector<double>& _inputPoints) {
        switch(m_functionType) {
            case FunctionType::SIN:
                return m_value * sin(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::COS:
                return m_value * cos(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::TAN:
                return m_value * tan(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::SINH:
                return m_value * sinh(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::COSH:
                return m_value * cosh(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::TANH:
                return m_value * tanh(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::LOG:
                return m_value * log(m_subfunctions[1]->evaluate(_inputPoints)) / log(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::LN:
                return m_value * log(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::POWER:
                return m_value * pow(m_subfunctions[0]->evaluate(_inputPoints), m_subfunctions[1]->evaluate(_inputPoints));
            case FunctionType::SQRT:
                return m_value * sqrt(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::EXP:
                return m_value * exp(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionType::CONSTANT:
                return m_value;
            case FunctionType::POWER_BASIS_POLYNOMIAL:
                break;
            case FunctionType::CHEBYSHEV_BASIS_POLYNOMIAL:
                break;
            case FunctionType::SUM:
                return sumEval(_inputPoints);
            case FunctionType::PRODUCT:
                return productEval(_inputPoints);
            case FunctionType::POWER_BASIS_MONOMIAL:
                return powerMonomialEval(_inputPoints);
            case FunctionType::CHEBYSHEV_BASIS_MONOMIAL:
                return chebyshevMonomialEval(_inputPoints);
            case FunctionType::CHEBYSHEV:
                return m_value * chebPower(m_subfunctions[0]->evaluate(_inputPoints), m_varIndex);
            case FunctionType::VARIABLE:
                return m_value * _inputPoints[m_varIndex];
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
        //Everything else must be a FunctionType string variable name, or subfunction name
                
        //Remove excess parenthesis surrounding the entire function.
        removeExtraParenthesis(_functionString);
        
        //TO start, split on +-.
        std::vector<std::string> subparts;
        splitSum(_functionString, subparts);
        
        if(subparts.size() == 0) {
            //This is an empty string! Something went wrong
            printAndThrowRuntimeError("Error Parsing Functions!");
        }
        else if(subparts.size() == 1) {
            //Reparse to look for product
            subparts.clear();
            splitProduct(_functionString, subparts, m_isMultiply);
            
            if(subparts.size() == 1) {
                //Parse the Coefficient
                parseCoeff(_functionString);

                //The whole thing is a more complex function
                parseComplexType(_functionString);
            }
            else {
                //Make this a PRODUCT type and split it up
                m_value = 1;
                parseProduct(subparts);
            }
        }
        else {
            //Make this a SUM type and split it up
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
                m_value += m_subfunctions[funcIndex]->getValue();
                m_subfunctions.erase(m_subfunctions.begin() + funcIndex);
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
                SIGNCHECKPRODUCT(m_isMultiply[funcIndex], m_value, m_subfunctions[funcIndex]->getValue());
                m_subfunctions.erase(m_subfunctions.begin() + funcIndex);
                m_isMultiply.erase(m_isMultiply.begin() + funcIndex);
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

    //Check for a constant LOG base
    if(m_functionType == FunctionType::LOG) {
        if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
            m_value /= log(m_subfunctions[0]->getValue());
            m_subfunctions.erase(m_subfunctions.begin());
            m_functionType = FunctionType::LN;
        }
    }
    
    //Check for a constant LN
    if(m_functionType == FunctionType::LN) {
        if(m_subfunctions[0]->getFunctionType() == FunctionType::CONSTANT) {
            m_value *= log(m_subfunctions[0]->getValue());
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
        
    //Check if a POWER is a POWER_BASIS_MONOMIAL
    /*if(m_functionType == FunctionType::POWER) {
        if(m_subfunctions[0]->getFunctionType() == FunctionType::VARIABLE && m_subfunctions[1]->getFunctionType() == FunctionType::CONSTANT) {
            size_t power = static_cast<size_t>(m_subfunctions[1]->getValue());
            if(power == m_subfunctions[0]->getValue()) {
                m_varIndexes.push_back(m_subfunctions[0]->getVarIndex());
                m_varPowers.push_back(power);
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionType::POWER_BASIS_MONOMIAL;
            }
        }
    }
    
    //Check if a CHEBYSHEV is a CHEBYSHEV_BASIS_MONOMIAL
    if(m_functionType == FunctionType::CHEBYSHEV) {
        if(m_subfunctions[0]->getFunctionType() == FunctionType::VARIABLE) {
            m_varIndexes.push_back(m_subfunctions[0]->getVarIndex());
            m_varPowers.push_back(m_varIndex);
            m_subfunctions.erase(m_subfunctions.begin());
            m_functionType = FunctionType::CHEBYSHEV_BASIS_MONOMIAL;
        }
    }*/
    
    //TODO: Fix this so we never change a child function. Ths causes issue as they are shared_ptrs and so will
    //be changing other functions as well.
    /*
    //Check if a PRODUCT has multiple POWER_BASIS_MONOMIAL or VARIABLE in it
    if(m_functionType == FunctionType::PRODUCT) {
        //Check if there is a first POWER_BASIS_MONOMIAL or VARIABLE.
        size_t funcIndex = 0;
        bool foundMonoimal = false;
        while(funcIndex < m_subfunctions.size()) {
            if((m_subfunctions[funcIndex]->getFunctionType() == FunctionType::POWER_BASIS_MONOMIAL || m_subfunctions[funcIndex]->getFunctionType() == FunctionType::VARIABLE) && m_isMultiply[funcIndex]) {
                foundMonoimal = true;
                break;
            }
            else {
                funcIndex++;
            }
        }
        
        //Merge other ones into it
        if(foundMonoimal) {
            size_t funcIndex2 = funcIndex + 1;
            while(funcIndex2 < m_subfunctions.size()) {
                if((m_subfunctions[funcIndex2]->getFunctionType() == FunctionType::POWER_BASIS_MONOMIAL || m_subfunctions[funcIndex2]->getFunctionType() == FunctionType::VARIABLE) && m_isMultiply[funcIndex2]) {
                    m_subfunctions[funcIndex]->mergePowerMonomialProduct(m_subfunctions[funcIndex2]);
                    m_subfunctions.erase(m_subfunctions.begin() + funcIndex2);
                }
                else {
                    funcIndex2++;
                }
            }
        }
        
        //If the size of the product is 1 and it had a POWER_BASIS_MONOMIAL it is now a POWER_BASIS_MONOMIAL
        if(m_subfunctions.size() == 1 && foundMonoimal) {
            m_value *= m_subfunctions[0]->getValue();
            m_varIndexes = m_subfunctions[0]->getVarIndexes();
            m_varPowers = m_subfunctions[0]->getVarPowers();
            m_subfunctions.erase(m_subfunctions.begin());
            m_functionType = FunctionType::POWER_BASIS_MONOMIAL;
        }
    }
    
    //Check if a PRODUCT has multiple CHEBYSHEV_BASIS_MONOMIAL in it
    if(m_functionType == FunctionType::PRODUCT) {
        //Check if there is a first CHEBYSHEV_BASIS_MONOMIAL.
        size_t funcIndex = 0;
        bool foundMonoimal = false;
        while(funcIndex < m_subfunctions.size()) {
            if(m_subfunctions[funcIndex]->getFunctionType() == FunctionType::CHEBYSHEV_BASIS_MONOMIAL && m_isMultiply[funcIndex]) {
                foundMonoimal = true;
                break;
            }
            else {
                funcIndex++;
            }
        }
        
        //Merge other ones into it
        if(foundMonoimal) {
            size_t funcIndex2 = funcIndex + 1;
            while(funcIndex2 < m_subfunctions.size()) {
                if(m_subfunctions[funcIndex2]->getFunctionType() == FunctionType::CHEBYSHEV_BASIS_MONOMIAL && m_isMultiply[funcIndex2]) {
                    //We can't merge CHEBYSHEV_BASIS_MONOMIAL if it has the same mons
                    if(m_subfunctions[funcIndex]->mergeChebyshevMonomialProduct(m_subfunctions[funcIndex2])) {
                        m_subfunctions.erase(m_subfunctions.begin() + funcIndex2);
                    }
                    else {
                        funcIndex2++;
                    }
                }
                else {
                    funcIndex2++;
                }
            }
        }
        
        //If the size of the product is 1 and it had a CHEBYSHEV_BASIS_MONOMIAL it is now a CHEBYSHEV_BASIS_MONOMIAL
        if(m_subfunctions.size() == 1 && foundMonoimal) {
            m_value *= m_subfunctions[0]->getValue();
            m_varIndexes = m_subfunctions[0]->getVarIndexes();
            m_varPowers = m_subfunctions[0]->getVarPowers();
            m_subfunctions.erase(m_subfunctions.begin());
            m_functionType = FunctionType::CHEBYSHEV_BASIS_MONOMIAL;
        }
    }
        */
    
    
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
    
    /*void mergePowerMonomialProduct(SharedFunctionPtr functionPtr) {
        //TODO: Delete this function?
        
        //Multiply two power basis monomials
        m_value *= functionPtr->getValue();
        std::vector<size_t>& otherPowers = functionPtr->getVarPowers();
        std::vector<size_t>& otherIndexes = functionPtr->getVarIndexes();
        
        for(size_t otherIndex = 0; otherIndex < otherPowers.size(); otherIndex++) {
            bool found = false;
            for(size_t thisIndex = 0; thisIndex < m_varPowers.size(); thisIndex++) {
                if(m_varIndexes[thisIndex] == otherIndexes[otherIndex]) {
                    m_varPowers[thisIndex] += otherPowers[otherIndex];
                    found = true;
                    break;
                }
            }
            if(!found) {
                m_varIndexes.push_back(otherIndexes[otherIndex]);
                m_varPowers.push_back(otherPowers[otherIndex]);
            }
        }
        //Turn a VARIABLE into a POWER_BASIS_MONOMIAL
        m_functionType = FunctionType::POWER_BASIS_MONOMIAL;
    }
    
    bool mergeChebyshevMonomialProduct(SharedFunctionPtr functionPtr) {
        //TODO: Delete this function?

        //Multiply two power basis monomials
        m_value *= functionPtr->getValue();
        std::vector<size_t>& otherPowers = functionPtr->getVarPowers();
        std::vector<size_t>& otherIndexes = functionPtr->getVarIndexes();
        
        std::vector<size_t> thisPowers = m_varPowers;
        std::vector<size_t> thisIndexes = m_varIndexes;

        for(size_t otherIndex = 0; otherIndex < otherPowers.size(); otherIndex++) {
            bool found = false;
            for(size_t thisIndex = 0; thisIndex < thisPowers.size(); thisIndex++) {
                if(thisIndexes[thisIndex] == otherIndexes[otherIndex]) {
                    return false;
                }
            }
            if(!found) {
                thisIndexes.push_back(otherIndexes[otherIndex]);
                thisPowers.push_back(otherPowers[otherIndex]);
            }
        }

        //Finalize
        m_varPowers = thisPowers;
        m_varIndexes = thisIndexes;
        return true;
    }
    */
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
    
    void splitSum(const std::string& _functionString, std::vector<std::string>& _subparts) {
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
            }
            else if(currChar == CHAR::LEFT_PAREN) {
                parenthesisCount++;
            }
            else if(parenthesisCount > 0) {
                //Do nothing, we are inside parenthesis
            }
            //Store the current substring if it exists
            else if(currChar == CHAR::PLUS || currChar == CHAR::MINUS) {
                if(currentSubstring.length() > 0) {
                    _subparts.push_back(currentSubstring);
                }
                currentSubstring = "";
            }
            currentSubstring += currChar;
        }
        
        //Add Anything Remaining
        _subparts.push_back(currentSubstring);
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
        //Add all the subfunctions
        for(size_t i = 0; i < _subparts.size(); i++) {
            addSubfunction(_subparts[i]);
        }
    }
    
    void parseCoeff(std::string& _functionString) {
        //Parses a number from the front of the string and removes it.
        std::string numStr = "";
        size_t toRemove = 0;
        size_t stringLength = _functionString.length();
        //Check if the next char is a number, - or .
        while(stringLength > toRemove) {
            if(isDigit(_functionString[toRemove])) {
                numStr += _functionString[toRemove++];
            }
            else if(_functionString[toRemove] == CHAR::PLUS) {
                toRemove++;
            }
            else {
                break;
            }
        }
        if(_functionString[toRemove] == CHAR::TIMES) {
            toRemove++;
        }
        
        //Get the truncated string
        _functionString = _functionString.substr(toRemove);
        //Parse the value.
        if(numStr == "") {
            m_value = 1.0;
        }
        else if(numStr == "-") {
            m_value = -1.0;
        }
        else {
            m_value = std::stod(numStr);
        }
    }
    
    void parseComplexType(std::string& _functionString) {
        size_t stringLenth = _functionString.length();
        
        //Parse Power
        if(parsePower(_functionString)) {
            m_functionType = FunctionType::POWER;
            return;
        }
        
        //Parse constants
        if(_functionString == "") {
            m_functionType = FunctionType::CONSTANT;
            return;
        }
        if(_functionString == "e" || _functionString == "E") {
            m_functionType = FunctionType::CONSTANT;
            m_value *= M_E;
            return;
        }
        if(_functionString == "pi" || _functionString == "PI" || _functionString == "Pi" || _functionString == "pI") {
            m_functionType = FunctionType::CONSTANT;
            m_value *= M_PI;
            return;
        }
        
        //Parse hyperbolic trig and SQRT
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
        }
        
        //Parse Normal Trig, EXP and LOG
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
                parseComplexFunctionParenthesis2(_functionString.substr(3));
                return;
            }
        }

        //Parse LN
        if(stringLenth > 2) {
            std::string lowercase2 = toLowerSubstring(_functionString, 0, 2);
            if(lowercase2 == "ln") {
                m_functionType = FunctionType::LN;
                parseComplexFunctionParenthesis(_functionString.substr(2));
                return;
            }
        }
        
        //Parse VARIABLE
        for(size_t i = 0; i < m_variableNames.size(); i++) {
            if(m_variableNames[i] == _functionString) {
                m_functionType = FunctionType::VARIABLE;
                m_varIndex = i;
                m_varIndexes.push_back(i);
                m_varPowers.push_back(1);
                return;
            }
        }
        
        //Parse CHEBYSHEV
        if(_functionString.length() > 1 && std::tolower(_functionString[0]) == 't' && isNumericDigit(_functionString[1])) {
            parseChebyshevFunction(_functionString);
            m_functionType = FunctionType::CHEBYSHEV;
            return;
        }
        
        printAndThrowRuntimeError("Function Type Unknown!");
    }
    
    void parseComplexFunctionParenthesis(const std::string& _functionString) {
        if(_functionString[0] != CHAR::LEFT_PAREN || _functionString[_functionString.length() - 1] != CHAR::RIGHT_PAREN) {
            printAndThrowRuntimeError("Failed to Parse Function!");
        }
        addSubfunction(_functionString.substr(1, _functionString.length() - 2));
    }

    void parseComplexFunctionParenthesis2(const std::string& _functionString) {
        if(_functionString[0] != CHAR::LEFT_PAREN || _functionString[_functionString.length() - 1] != CHAR::RIGHT_PAREN) {
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
            else if(c == CHAR::LEFT_PAREN) {
                parenthesisCount++;
            }
            //Check for a comma and add to the right substring
            if(c == CHAR::COMMA && parenthesisCount == 0) {
                part2 = true;
            }
            else if (part2) {
                substring2 += c;
            }
            else {
                substring1 += c;
            }
        }
        addSubfunction(substring1);
        addSubfunction(substring2);
    }

    void parseChebyshevFunction(const std::string& _functionString) {
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

    bool parsePower(const std::string& _functionString) {
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
            //Check for a comma and add to the right substring
            if(c == CHAR::POWER && parenthesisCount == 0) {
                part2 = true;
            }
            else if (part2) {
                substring2 += c;
            }
            else {
                substring1 += c;
            }
        }
        if(part2) {
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
    double powerMonomialEval(const std::vector<double>& _inputPoints) {
        double result = m_value;
        for(size_t i = 0; i < m_varIndexes.size(); i++) {
            result *= power(_inputPoints[m_varIndexes[i]], m_varPowers[i]);
        }
        return result;
    }

    double chebyshevMonomialEval(const std::vector<double>& _inputPoints) {
        double result = m_value;
        for(size_t i = 0; i < m_varIndexes.size(); i++) {
            result *= chebPower(_inputPoints[m_varIndexes[i]], m_varPowers[i]);
        }
        return result;
    }
    
    double sumEval(const std::vector<double>& _inputPoints) {
        double result = m_value;
        for(size_t i = 0; i < m_subfunctions.size(); i++) {
            result += m_subfunctions[i]->evaluate(_inputPoints);
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
    
    double productEval(const std::vector<double>& _inputPoints) {
        double result = m_value * m_subfunctions[0]->evaluate(_inputPoints);
        for(size_t i = 1; i < m_isMultiply.size(); i++) {
            SIGNCHECKPRODUCT(m_isMultiply[i], result, m_subfunctions[i]->evaluate(_inputPoints));
        }
        return result;
    }
    
public:
    void getFunctionLevels(std::vector<std::vector<SharedFunctionPtr>>& allFunctionLevels) {
        //The function level is 1 + max(function level of subfunctions).
        //If a function has no subfunctions, the function level is 0.
        //This way when we evaluate we can evaluate all the functions starting a level 0 at going up.
        
        m_functionLevel = 0;
        for(size_t i = 0; i < m_subfunctions.size(); i++) {
            m_subfunctions[i]->getFunctionLevels(allFunctionLevels);
            size_t currLevel = m_subfunctions[i]->getFunctionLevel();
            m_functionLevel = std::max(m_functionLevel, currLevel + 1);
            
            if(allFunctionLevels.size() <= currLevel) {
                allFunctionLevels.resize(currLevel + 1);
            }
            allFunctionLevels[currLevel].push_back(m_subfunctions[i]);
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
    std::vector<size_t>& getVarIndexes() {
        return m_varIndexes;
    }
    
    std::vector<size_t>& getVarPowers() {
        return m_varPowers;
    }

    const std::vector<bool>& getIsMultiply() const {
        return m_isMultiply;
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
    
    const std::vector<size_t>& getVarIndexes() const {
        return m_varIndexes;
    }
    
    const std::vector<size_t>& getVarPowers() const {
        return m_varPowers;
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

    const size_t getDimension() const {
        return m_dimension;
    }

    const size_t getNumUsedDimensions() const {
        return m_numUsedDimensions;
    }

    const std::vector<EvaluateGridInfo>& getEvaluateGridInfo() const {
        return m_evaluteGridInfo;
    }

    const EvaluateGridType getEvaluateGridType() const {
        return m_evaluateGridType;
    }

    const std::vector<std::vector<SharedFunctionPtr>>& getAllFunctionLevels() const {
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
                SharedFunctionPtr currFunction = std::make_shared<Function>(*it->second.get(), true);
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
    Function(const Function& other, bool isShallowCopy = false) :
    m_isTopFunction(other.isTopFunction()),
    m_functionLevel(other.getFunctionLevel()),
    m_isMultiply(other.getIsMultiply()),
    m_functionType(other.getFunctionType()),
    m_value(other.getValue()),
    m_varIndex(other.getVarIndex()),
    m_varIndexes(other.getVarIndexes()),
    m_varPowers(other.getVarPowers()),
    m_functionString(other.getFunctionString()),
    m_functionName(other.getFunctionName()),
    m_variableNames(other.getVariableNames()),
    m_dimension(other.getDimension()),
    m_hasDimension(other.getHasDimension()),
    m_numUsedDimensions(other.getNumUsedDimensions()),
    m_partialEvaluations(other.getPartialEvals()),
    m_evaluteGridInfo(other.getEvaluateGridInfo()),
    m_evaluateGridType(other.getEvaluateGridType())
    {
        //Copies everything except the pointers
        if(!isShallowCopy) {
            printAndThrowRuntimeError("Can only do a shallow copy of functions!");
        }
    }
    
    //Helper to the copy constructor, updates the pointers
    void finalizeCopy(const Function& other) {
        const size_t threadNum = s_allFunctions.size() - 1; //Copy from this spot in the s_allFunctions
        
        //Copy m_allFunctionLevels
        const std::vector<std::vector<SharedFunctionPtr>>& functionLevelsToCopy = other.getAllFunctionLevels();
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
    
public:
    //A function should only be created through this. It adds them to the maps correctly.
    static SharedFunctionPtr addFunction(const std::string& _functionName, const std::string& _functionString, const std::vector<std::string>& _variableNames) {
        if(unlikely(s_allFunctions.size() == 0)) {
            s_allFunctions.resize(1);
            s_allFunctionNames.resize(1);
        }

        SharedFunctionPtr function = std::make_shared<Function>(_functionName, _functionString, _variableNames);
        
        //Add to the maps
        if(s_allFunctions[0].find(_functionString) != s_allFunctions[0].end()) {
            printAndThrowRuntimeError("Trying to add function Twice! Name is " + _functionName + " String is " + _functionString);
        }
        else {
            s_allFunctions[0][_functionString] = function;
        }
        
        if(_functionName != "") { //Only track it if it has a name
            if(s_allFunctionNames[0].find(_functionName) != s_allFunctionNames[0].end()) {
                printAndThrowRuntimeError("Trying to add function Twice! Name is " + _functionName + " String is " + _functionString);
            }
            else {
                s_allFunctionNames[0][_functionName] = function;
            }
        }
        return function;
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
    std::vector<std::vector<SharedFunctionPtr>> m_allFunctionLevels;

    //The subfunctions and coresponding signs
    std::vector<SharedFunctionPtr>              m_subfunctions;
    std::vector<bool>                           m_isMultiply;
    
    //The type of function
    FunctionType                               m_functionType;
    
    //The Signed Coefficient of anything that isn't a sum. Sums all have coefficient of 1.
    double                                      m_value;
    
    //For VARIABLE: The variable index
    size_t                                      m_varIndex;
    
    //For POWER_BASIS_MONOMIAL and CHEBYSHEV_BASIS_MONOMIAL: The variable indexes
    std::vector<size_t>                         m_varIndexes;
    std::vector<size_t>                         m_varPowers;
        
    //Variable Names and Subfunctions Names
    std::string                                 m_functionString;
    std::string                                 m_functionName;
    std::vector<std::string>                    m_variableNames;
    
    //For evaluate grid
    size_t                                      m_dimension;
    std::vector<bool>                           m_hasDimension;
    size_t                                      m_numUsedDimensions;
    std::vector<double>                         m_partialEvaluations;
    std::vector<EvaluateGridInfo>               m_evaluteGridInfo;
    EvaluateGridType                            m_evaluateGridType;
    
    static size_t           m_timerEvaluateGridIndex;
    Timer&                  m_timer = Timer::getInstance();
};

std::vector<Function::FunctionMap>    Function::s_allFunctions;
std::vector<Function::FunctionMap>    Function::s_allFunctionNames;
size_t                      Function::m_timerEvaluateGridIndex = -1;


#endif /* Function_h */
