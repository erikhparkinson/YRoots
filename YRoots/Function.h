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

//TODO: Make it so any complex function can handle of coefficient
//Sums can have coefficients in from of them as well.
//That way we can remove the signs, it will just be coefficients instead
//Products Still have to have the sign variable to know if it will divide or multiply

//When Parsing a complex thing, first go as far as possible while it is numbers or '.'
//Then check the type.

enum FunctionTypes{
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
    SUBFUNCTION
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

#define SIGNCHECKPRODUCT(isPositive, sum, number) (isPositive ? sum *= number : sum /= number)

//TODO: This class should really be function part.
//Have Another class the holds the subfunctions of the functions and the function itself.
    //It has to figure out which of the subfunctions are actually in the function
//When Evaluate is Called, it first evaluates all of the subfunctions. Then it evaluates the actual function.
//The subfunction should have a seperate pre_evaluate function call that actually computes things, and then
//have evaluate just return the result.

class Function{
public:
    Function(std::string _functionString, std::vector<std::string> _variableNames, std::vector<std::string> _subfunctionNames) :
    m_functionType(FunctionTypes::CONSTANT),
    m_value(0),
    m_varIndex(0),
    m_subfunctionEval(0),
    m_isSubfunction(false)
    {
        m_functionString = _functionString;
        m_variableNames = _variableNames;
        m_subfunctionNames = _subfunctionNames;
        functionParse(_functionString, _variableNames, _subfunctionNames);
    }
    
    Function(const Function& other)
    {
        copyFunction(other);
    }
        
    ~Function(){
        
    }
        
    double preEvaluate(const std::vector<double>& _inputPoints) {
        assert(m_isSubfunction);
        m_isSubfunction = false;
        double returnVal = evaluate(_inputPoints);
        m_isSubfunction = true;
        return returnVal;
    }

    void evaluateGrid(const std::vector<std::vector<double>>& _grid, double* _results, double _divisor = 1.0)
    {
        //Every Function will have to have a vector created inside of a certain size based on the max size of the grid.
        //Evaluate Grid will then pass evaluate grid into it's subfunctions if it has any. Polynomials, constants, and monomials
        //Will just evaluate the grid and store them in the vector. Then higher level functions iterate throug these vectors and
        //run their function on them.
        
        //So really double* _results will only get called on the top function class, every other class just gets _grid
        //passed into it and evaluates into the vector. Then the top function calss writes it's result into _results.
        
        //And grid is only used in the bottom level functions, everything else just uses the vectors in the functions below it.
        //These are now non-grid like so all the points just have to be evaluated.
        
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
        _results[evalSpot++] = evaluate(inputPoints)/_divisor;
        while (spotToInc < dimension) {
            bool firstPass = true;
            while(++inputSpot[spotToInc] < numPoints) {
                inputPoints[spotToInc] = _grid[spotToInc][inputSpot[spotToInc]];
                _results[evalSpot++] = evaluate(inputPoints)/_divisor;
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
        if(unlikely(m_isSubfunction)) {
            return m_subfunctionEval;
        }

        switch(m_functionType) {
            case FunctionTypes::SIN:
                return m_value * sin(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::COS:
                return m_value * cos(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::TAN:
                return m_value * tan(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::SINH:
                return m_value * sinh(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::COSH:
                return m_value * cosh(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::TANH:
                return m_value * tanh(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::LOG:
                return m_value * log(m_subfunctions[1]->evaluate(_inputPoints)) / log(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::LN:
                return m_value * log(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::POWER:
                return m_value * pow(m_subfunctions[0]->evaluate(_inputPoints), m_subfunctions[1]->evaluate(_inputPoints));
            case FunctionTypes::SQRT:
                return m_value * sqrt(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::EXP:
                return m_value * exp(m_subfunctions[0]->evaluate(_inputPoints));
            case FunctionTypes::CONSTANT:
                return m_value;
            case FunctionTypes::POWER_BASIS_POLYNOMIAL:
                break;
            case FunctionTypes::CHEBYSHEV_BASIS_POLYNOMIAL:
                break;
            case FunctionTypes::SUM:
                return sumEval(_inputPoints);
            case FunctionTypes::PRODUCT:
                return productEval(_inputPoints);
            case FunctionTypes::POWER_BASIS_MONOMIAL:
                return powerMonomialEval(_inputPoints);
            case FunctionTypes::CHEBYSHEV_BASIS_MONOMIAL:
                return chebyshevMonomialEval(_inputPoints);
            case FunctionTypes::CHEBYSHEV:
                return m_value * chebPower(m_subfunctions[0]->evaluate(_inputPoints), m_varIndex);
            case FunctionTypes::VARIABLE:
                return m_value * _inputPoints[m_varIndex];
            default:
                throw std::runtime_error("Unknown Function Type Encountered in evaluate! Fix Switch Statement!");
                break;
        }
        //TODO: Get rid of this Once everything else is implemented.
        return 0.0;
    }
        
private:
    void functionParse(std::string& _functionString, std::vector<std::string>& _variableNames, std::vector<std::string>& _subfunctionNames) {
        //GRAMMER:
        // ( and ) make subfunctions
        // +,-,*,/,^ are math symbols
        //numbers, e, pi are constants
        //Everything else must be a FunctionTypes string variable name, or subfunction name
                
        //Remove excess parenthesis surrounding the entire function.
        removeExtraParenthesis(_functionString);
        
        //TO start, split on +-.
        std::vector<std::string> subparts;
        splitSum(_functionString, subparts);
        
        if(subparts.size() == 0) {
            //This is an empty string! Something went wrong
            std::string errorMessage = "Error Parsing Functions!";
            std::cout<<errorMessage<<"\n";
            throw std::runtime_error(errorMessage);
        }
        else if(subparts.size() == 1) {
            //Parse the Coefficient
            parseCoeff(_functionString);
            
            //Reparse to look for product
            subparts.clear();
            splitProduct(_functionString, subparts, m_isMultiply);
            
            if(subparts.size() == 1) {
                //The whole thing is a more complex function
                parseComplexType(_functionString, _variableNames, _subfunctionNames);
            }
            else {
                //Make this a PRODUCT type and split it up
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
    
    void simplifyExpressions()
    {
        //Pull any constants out of a sum
        if(m_functionType == FunctionTypes::SUM) {
            //Pull any constants out of a sum
            size_t funcIndex = 0;
            while(funcIndex < m_subfunctions.size()) {
                if(m_subfunctions[funcIndex]->getFunctionType() == FunctionTypes::CONSTANT) {
                    m_value += m_subfunctions[funcIndex]->getValue();
                    m_subfunctions.erase(m_subfunctions.begin() + funcIndex);
                }
                else {
                    funcIndex++;
                }
            }
            //If the size of the sum is 0 it is now a constant
            if(m_subfunctions.size() == 0) {
                m_functionType = FunctionTypes::CONSTANT;
            }
        }
        
        //Pull any constants out of a product. Check if we multiply or divide!
        if(m_functionType == FunctionTypes::PRODUCT) {
            //Pull any constants out of a product
            size_t funcIndex = 0;
            while(funcIndex < m_subfunctions.size()) {
                if(m_subfunctions[funcIndex]->getFunctionType() == FunctionTypes::CONSTANT) {
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
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a constant LOG base
        if(m_functionType == FunctionTypes::LOG) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value /= log(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::LN;
            }
        }
        
        //Check for a constant LN
        if(m_functionType == FunctionTypes::LN) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= log(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant SIN
        if(m_functionType == FunctionTypes::SIN) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= sin(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant COS
        if(m_functionType == FunctionTypes::COS) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= cos(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant TAN
        if(m_functionType == FunctionTypes::TAN) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= tan(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant SINH
        if(m_functionType == FunctionTypes::SINH) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= sinh(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant COSH
        if(m_functionType == FunctionTypes::COSH) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= cosh(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant TANH
        if(m_functionType == FunctionTypes::TANH) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= tanh(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant SQRT
        if(m_functionType == FunctionTypes::SQRT) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= sqrt(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant EXP
        if(m_functionType == FunctionTypes::EXP) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= exp(m_subfunctions[0]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }

        //Check for a contant POWER
        if(m_functionType == FunctionTypes::POWER) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT && m_subfunctions[1]->getFunctionType() == FunctionTypes::CONSTANT) {
                m_value *= pow(m_subfunctions[0]->getValue(), m_subfunctions[1]->getValue());
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
            }
        }
        
        //Check powers with ones and zeros
        if(m_functionType == FunctionTypes::POWER) {
            //Check if it is something to the POWER of 0
            if(m_subfunctions[1]->getFunctionType() == FunctionTypes::CONSTANT && m_subfunctions[1]->getValue() == 0.0) {
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
                m_value = 1;
            }
            //Check if it is something to the POWER of 1
            else if(m_subfunctions[1]->getFunctionType() == FunctionTypes::CONSTANT && m_subfunctions[1]->getValue() == 1.0) {
                std::shared_ptr tempFunctionPtr(m_subfunctions[0]);
                copyFunction(*tempFunctionPtr.get());
            }
            //Check it is 0 to the POWER of something
            else if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT && m_subfunctions[0]->getValue() == 0.0) {
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
                m_value = 0;
            }
            //Check it is 1 to the POWER of something
            else if(m_subfunctions[0]->getFunctionType() == FunctionTypes::CONSTANT && m_subfunctions[0]->getValue() == 1.0) {
                m_subfunctions.erase(m_subfunctions.begin());
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CONSTANT;
                m_value = 1;
            }
        }

        //Check if a POWER is a POWER_BASIS_MONOMIAL
        if(m_functionType == FunctionTypes::POWER) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::VARIABLE && m_subfunctions[1]->getFunctionType() == FunctionTypes::CONSTANT) {
                size_t power = static_cast<size_t>(m_subfunctions[1]->getValue());
                if(power == m_subfunctions[0]->getValue()) {
                    m_varIndexes.push_back(m_subfunctions[0]->getVarIndex());
                    m_varPowers.push_back(power);
                    m_subfunctions.erase(m_subfunctions.begin());
                    m_subfunctions.erase(m_subfunctions.begin());
                    m_functionType = FunctionTypes::POWER_BASIS_MONOMIAL;
                }
            }
        }
        
        //Check if a CHEBYSHEV is a CHEBYSHEV_BASIS_MONOMIAL
        if(m_functionType == FunctionTypes::CHEBYSHEV) {
            if(m_subfunctions[0]->getFunctionType() == FunctionTypes::VARIABLE) {
                m_varIndexes.push_back(m_subfunctions[0]->getVarIndex());
                m_varPowers.push_back(m_varIndex);
                m_subfunctions.erase(m_subfunctions.begin());
                m_functionType = FunctionTypes::CHEBYSHEV_BASIS_MONOMIAL;
            }
        }
        
        //Check if a PRODUCT has multiple POWER_BASIS_MONOMIAL or VARIABLE in it
        if(m_functionType == FunctionTypes::PRODUCT) {
            //Check if there is a first POWER_BASIS_MONOMIAL or VARIABLE.
            size_t funcIndex = 0;
            bool foundMonoimal = false;
            while(funcIndex < m_subfunctions.size()) {
                if((m_subfunctions[funcIndex]->getFunctionType() == FunctionTypes::POWER_BASIS_MONOMIAL || m_subfunctions[funcIndex]->getFunctionType() == FunctionTypes::VARIABLE) && m_isMultiply[funcIndex]) {
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
                    if((m_subfunctions[funcIndex2]->getFunctionType() == FunctionTypes::POWER_BASIS_MONOMIAL || m_subfunctions[funcIndex2]->getFunctionType() == FunctionTypes::VARIABLE) && m_isMultiply[funcIndex2]) {
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
                m_functionType = FunctionTypes::POWER_BASIS_MONOMIAL;
            }
        }
        
        //Check if a PRODUCT has multiple CHEBYSHEV_BASIS_MONOMIAL in it
        if(m_functionType == FunctionTypes::PRODUCT) {
            //Check if there is a first CHEBYSHEV_BASIS_MONOMIAL.
            size_t funcIndex = 0;
            bool foundMonoimal = false;
            while(funcIndex < m_subfunctions.size()) {
                if(m_subfunctions[funcIndex]->getFunctionType() == FunctionTypes::CHEBYSHEV_BASIS_MONOMIAL && m_isMultiply[funcIndex]) {
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
                    if(m_subfunctions[funcIndex2]->getFunctionType() == FunctionTypes::CHEBYSHEV_BASIS_MONOMIAL && m_isMultiply[funcIndex2]) {
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
                m_functionType = FunctionTypes::CHEBYSHEV_BASIS_MONOMIAL;
            }
        }

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
        

        //TODO: Check the subfunctions are SUBFUNCTION, in that case replace them with the correct shared_ptr
        
        //TODO: At the end go through and see if any Functions are equal to each other, if so make them the same shared ptr.
        //Make a Tree, constants and variables are level 0 and everything else is level 1 + max of sublevels.
        //When evaluating go through and evaluate each function in level 0 then 1, etc.
        //For evaluate grid pass in the size of the max vector needed. Then have each function store the values in a grid.
        
        //TODO: There should be no constants left at the end of this function! That should all be assimilated into the others.
        //Unless the whole function is one big constant I guess.
    }
    
    void mergePowerMonomialProduct(std::shared_ptr<Function> functionPtr) {
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
        m_functionType = FunctionTypes::POWER_BASIS_MONOMIAL;
    }
    
    bool mergeChebyshevMonomialProduct(std::shared_ptr<Function> functionPtr) {
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
                if(currentSubstring.length() > 0) {
                    _subparts.push_back(currentSubstring);
                    _isMultiply.push_back(currentIsMultiply);
                }
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
        m_functionType = FunctionTypes::SUM;
        m_value = 0;
        //Add all the subfunctions
        for(size_t i = 0; i < _subparts.size(); i++) {
            m_subfunctions.push_back(std::make_shared<Function>(_subparts[i], m_variableNames, m_subfunctionNames));
        }
    }

    void parseProduct(const std::vector<std::string>& _subparts) {
        m_functionType = FunctionTypes::PRODUCT;
        //Add all the subfunctions
        for(size_t i = 0; i < _subparts.size(); i++) {
            m_subfunctions.push_back(std::make_shared<Function>(_subparts[i], m_variableNames, m_subfunctionNames));
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
    
    void parseComplexType(std::string& _functionString, std::vector<std::string>& _variableNames, std::vector<std::string>& _subfunctionNames) {
        size_t stringLenth = _functionString.length();
        
        //Parse Power
        if(parsePower(_functionString)) {
            m_functionType = FunctionTypes::POWER;
            return;
        }
        
        //Parse constants
        if(_functionString == "") {
            m_functionType = FunctionTypes::CONSTANT;
            return;
        }
        if(_functionString == "e" || _functionString == "E") {
            m_functionType = FunctionTypes::CONSTANT;
            m_value *= M_E;
            return;
        }
        if(_functionString == "pi" || _functionString == "PI" || _functionString == "Pi" || _functionString == "pI") {
            m_functionType = FunctionTypes::CONSTANT;
            m_value *= M_PI;
            return;
        }
        
        //Parse hyperbolic trig and SQRT
        if(stringLenth > 4) {
            std::string lowercase4 = toLowerSubstring(_functionString, 0, 4);
            if(lowercase4 == "sinh") {
                m_functionType = FunctionTypes::SINH;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "cosh") {
                m_functionType = FunctionTypes::COSH;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "tanh") {
                m_functionType = FunctionTypes::TANH;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
            if(lowercase4 == "sqrt") {
                m_functionType = FunctionTypes::SQRT;
                parseComplexFunctionParenthesis(_functionString.substr(4));
                return;
            }
        }
        
        //Parse Normal Trig, EXP and LOG
        if(stringLenth > 3) {
            std::string lowercase3 = toLowerSubstring(_functionString, 0, 3);
            if(lowercase3 == "sin") {
                m_functionType = FunctionTypes::SIN;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "cos") {
                m_functionType = FunctionTypes::COS;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "tan") {
                m_functionType = FunctionTypes::TAN;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "exp") {
                m_functionType = FunctionTypes::EXP;
                parseComplexFunctionParenthesis(_functionString.substr(3));
                return;
            }
            if(lowercase3 == "log") {
                m_functionType = FunctionTypes::LOG;
                parseComplexFunctionParenthesis2(_functionString.substr(3));
                return;
            }
        }

        //Parse LN
        if(stringLenth > 2) {
            std::string lowercase2 = toLowerSubstring(_functionString, 0, 2);
            if(lowercase2 == "ln") {
                m_functionType = FunctionTypes::LN;
                parseComplexFunctionParenthesis(_functionString.substr(2));
                return;
            }
        }
        
        //Parse VARIABLE
        for(size_t i = 0; i < _variableNames.size(); i++) {
            if(_variableNames[i] == _functionString) {
                m_functionType = FunctionTypes::VARIABLE;
                m_varIndex = i;
                m_varIndexes.push_back(i);
                m_varPowers.push_back(1);
                return;
            }
        }
        
        //Parse SUBFUNCTIONS
        for(size_t i = 0; i < _subfunctionNames.size(); i++) {
            if(_subfunctionNames[i] == _functionString) {
                m_functionType = FunctionTypes::SUBFUNCTION;
                m_varIndex = i;
                return;
            }
        }

        //Parse CHEBYSHEV
        if(_functionString.length() > 1 && std::tolower(_functionString[0]) == 't' && isNumericDigit(_functionString[1])) {
            parseChebyshevFunction(_functionString);
            m_functionType = FunctionTypes::CHEBYSHEV;
            return;
        }
        
        printAndThrowRuntimeError("Function Type Unknown!");
    }
    
    void parseComplexFunctionParenthesis(const std::string& _functionString) {
        if(_functionString[0] != CHAR::LEFT_PAREN || _functionString[_functionString.length() - 1] != CHAR::RIGHT_PAREN) {
            printAndThrowRuntimeError("Failed to Parse Function!");
        }
        m_subfunctions.push_back(std::make_shared<Function>(_functionString.substr(1, _functionString.length() - 2), m_variableNames, m_subfunctionNames));
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
        
        m_subfunctions.push_back(std::make_shared<Function>(substring1, m_variableNames, m_subfunctionNames));
        m_subfunctions.push_back(std::make_shared<Function>(substring2, m_variableNames, m_subfunctionNames));
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
        m_subfunctions.push_back(std::make_shared<Function>(substring2, m_variableNames, m_subfunctionNames));
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
            m_subfunctions.push_back(std::make_shared<Function>(substring1, m_variableNames, m_subfunctionNames));
            m_subfunctions.push_back(std::make_shared<Function>(substring2, m_variableNames, m_subfunctionNames));
            return true;
        }
        return false;
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
    
//Getters
public:
    std::vector<size_t>& getVarIndexes() {
        return m_varIndexes;
    }
    
    std::vector<size_t>& getVarPowers() {
        return m_varPowers;
    }

    std::vector<std::shared_ptr<Function>> getSubfunctions() const {
        return m_subfunctions;
    }

    std::vector<bool> getIsMultiply() const {
        return m_isMultiply;
    }

    FunctionTypes getFunctionType() const {
        return m_functionType;
    }
    
    double getValue() const {
        return m_value;
    }
    
    size_t getVarIndex() const {
        return m_varIndex;
    }
    
    std::vector<size_t> getVarIndexes() const {
        return m_varIndexes;
    }
    
    std::vector<size_t> getVarPowers() const {
        return m_varPowers;
    }

    double getSubfunctionEval() const {
        return m_subfunctionEval;
    }
    
    bool getIsSubfunction() const {
        return m_isSubfunction;
    }

    std::string getFunctionString() const {
        return m_functionString;
    }

    std::vector<std::string> getVariableNames() const {
        return m_variableNames;
    }

    std::vector<std::string> getSubfunctionNames() const {
        return m_subfunctionNames;
    }
    
private:
    void copyFunction(const Function& other)
    {
        m_subfunctions = other.getSubfunctions();
        m_isMultiply = other.getIsMultiply();
        m_functionType = other.getFunctionType();
        m_value = other.getValue();
        m_varIndex = other.getVarIndex();
        m_varIndexes = other.getVarIndexes();
        m_varPowers = other.getVarPowers();
        m_subfunctionEval = other.getSubfunctionEval();
        m_isSubfunction = other.getIsSubfunction();
        m_functionString = other.getFunctionString();
        m_variableNames = other.getVariableNames();
        m_subfunctionNames = other.getSubfunctionNames();
    }
    
private:
    //The subfunctions and coresponding signs
    std::vector<std::shared_ptr<Function>>      m_subfunctions;
    std::vector<bool>                           m_isMultiply;
    
    //The type of function
    FunctionTypes                               m_functionType;
    
    //The Signed Coefficient of anything that isn't a sum. Sums all have coefficient of 1.
    double                                      m_value;
    
    //For VARIABLE: The variable index
    size_t                                      m_varIndex;
    
    //For POWER_BASIS_MONOMIAL and CHEBYSHEV_BASIS_MONOMIAL: The variable indexes
    std::vector<size_t>                         m_varIndexes;
    std::vector<size_t>                         m_varPowers;
    
    //For Subfunctions
    double                                      m_subfunctionEval;
    bool                                        m_isSubfunction;
    
    //Variable Names and Subfunctions Names
    std::string                                 m_functionString;
    std::vector<std::string>                    m_variableNames;
    std::vector<std::string>                    m_subfunctionNames;
};



























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
    
    virtual void evaluateGrid(const std::vector<std::vector<double>>& grid, double* results, double divisor = 1.0)
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
        results[evalSpot++] = evaluate(inputPoints)/divisor;
        while (spotToInc < dimension) {
            bool firstPass = true;
            while(++inputSpot[spotToInc] < numPoints) {
                inputPoints[spotToInc] = grid[spotToInc][inputSpot[spotToInc]];
                results[evalSpot++] = evaluate(inputPoints)/divisor;
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
