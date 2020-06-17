//
//  PowerBasisPolynomial.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef PowerBasisPolynomial_h
#define PowerBasisPolynomial_h

#include "Function.h"
#include <iostream>

class PowerBasisPolynomial: public FunctionInterface {
public:
    PowerBasisPolynomial(const std::string& functionString, const std::vector<std::string>& variableNames) :
    FunctionInterface(functionString, variableNames)
    {
        m_dim = variableNames.size();
        
        std::vector<std::pair<double, std::vector<size_t>>> coeffs;
        
        //Get the monmials
        size_t last = 0;
        size_t next = 0;
        double sign = 1;
        while (true) {
            size_t nextPlus = functionString.find("+", last);
            size_t nextMinus = functionString.find("-", last);
            next = std::min(nextPlus, nextMinus);
            if(next == std::string::npos) {
                break;
            }
            
            //Get the coefficient
            std::string monomial = functionString.substr(last, next-last);
            if(monomial!="") {
                coeffs.push_back(getMonomial(sign, monomial, variableNames));
            }
            
            last = next + 1;
            //Get the sign
            if(nextPlus < nextMinus) {
                sign = 1;
            }
            else {
                sign = -1;
            }
        }
        //Get the coefficient
        std::string monomial = functionString.substr(last, next-last);
        coeffs.push_back(getMonomial(sign, monomial, variableNames));
        
        //Get the dimensions
        for(size_t i = 0; i < m_dim; i++) {
            m_matrixDim.push_back(0);
            for(size_t j = 0; j < coeffs.size(); j++) {
                m_matrixDim[i] = std::max(m_matrixDim[i], coeffs[j].second[i] + 1);
            }
        }
        
        //Remake the m_matrixDim
        std::vector<size_t> oldMatrixDim = m_matrixDim;
        m_matrixDim.clear();
        for(size_t i = 0; i < m_dim; i++) {
            m_matrixDim.push_back(oldMatrixDim[m_dim - i - 1]);
        }
        
        //Get the matrixDimProdects
        for(size_t i = 0; i < m_dim; i++) {
            size_t newProduct = 1;
            for(size_t j = 0; j < i; j++) {
                newProduct *= m_matrixDim[j];
            }
            m_matrixDimProducts.push_back(newProduct);
        }
        
        //Put everything into an array
        size_t arraySize = 1;
        for(size_t i = 0; i < m_matrixDim.size(); i++) {
            arraySize *= m_matrixDim[i];
        }
        
        //Allocate the Array Memory
        double* coeffArray = (double*) malloc(arraySize * sizeof (double));
        if (!coeffArray) {
            throw std::runtime_error("Memory Allocation Failed!");
        }
        
        //Zero out the array
        for(size_t i = 0; i < arraySize; i++) {
            coeffArray[i] = 0;
        }
        
        //Put the data in the array
        for(size_t i = 0; i < coeffs.size(); i++) {
            coeffArray[getArraySpot(coeffs[i].second)] = coeffs[i].first;
        }
                
        m_arraySizes.push_back(arraySize);
        m_arrays.push_back(coeffArray);
        
        //Create the other arrays needed for quick evaluation
        size_t reductionSize = arraySize;
        for(size_t i = 0; i < m_dim; i++) {
            reductionSize /= m_matrixDim[i];
            m_arraySizes.push_back(reductionSize);
            double* reductionArray = (double*) malloc(reductionSize * sizeof (double));
            if (!reductionArray) {
                free(reductionArray);
                throw std::runtime_error("Memory Allocation Failed!");
            }
            m_arrays.push_back(reductionArray);
        }
    }
    
    PowerBasisPolynomial(PowerBasisPolynomial const&) = delete;
    PowerBasisPolynomial& operator=(PowerBasisPolynomial const&) = delete;
    
    ~PowerBasisPolynomial()
    {
        for(size_t i = 0; i < m_arrays.size(); i++){
            free(m_arrays[i]);
        }
    }
    
    std::pair<double, std::vector<size_t>> getMonomial(double sign, std::string monomial, const std::vector<std::string>& variableNames) {
        double coeff = 1;
        std::vector<std::string> variablesParts = split(monomial, "*");
        if(variablesParts.size() == 0) {
            throw std::runtime_error("Failure parsing function!");
        }
        
        std::vector<std::vector<std::string>> monomialPowers;
        for(size_t coeffIndex = 0; coeffIndex < variablesParts.size(); coeffIndex++) {
            if(is_number(variablesParts[coeffIndex])) {
                coeff = std::stod(variablesParts[coeffIndex]);
            }
            else {
                monomialPowers.push_back(split(variablesParts[coeffIndex], "^"));
            }
        }
        
        //Get the powers
        std::vector<size_t> powers(m_dim, 0);
        for(size_t variableIndex = 0; variableIndex < m_dim; variableIndex++) {
            std::string currVariableName = variableNames[variableIndex];
            for(size_t powerIndex = 0; powerIndex < monomialPowers.size(); powerIndex++) {
                if(monomialPowers[powerIndex][0] == currVariableName) {
                    if(monomialPowers[powerIndex].size() > 1) {
                        powers[variableIndex] = std::stod(monomialPowers[powerIndex][1]);
                    }
                    else {
                        powers[variableIndex] = 1.0;
                    }
                    break;
                }
            }
        }
        
        return std::make_pair(coeff*sign, powers);
    }
    
    size_t getArraySpot(std::vector<size_t> spots) {
        size_t val = 0;
        for(size_t i = 0; i < spots.size(); i++) {
            val += spots[m_dim - i - 1] * m_matrixDimProducts[i];
        }
        return val;
    }
    
    double evaluate(const std::vector<double>& inputPoints) override {
        size_t currDim = 0;
        while(currDim < m_dim) {
            m_currPoint = inputPoints[m_dim - currDim - 1];
            size_t currArraySpot = m_arraySizes[currDim] - 1;
            int64_t nextArraySpot = m_arraySizes[currDim+1] - 1;
            while(nextArraySpot >= 0) {
                m_currDouble = m_arrays[currDim][currArraySpot];
                for(size_t i = 0; i + 1 < m_matrixDim[currDim]; i++) {
                    currArraySpot--;
                    m_currDouble = (m_currDouble*m_currPoint) + m_arrays[currDim][currArraySpot];
                }
                m_arrays[currDim + 1][nextArraySpot] = m_currDouble;
                currArraySpot--;
                nextArraySpot--;
            }
            currDim++;
        }
        return m_arrays[m_dim][0];
    }
    
private:
    size_t                  m_dim;
    std::vector<size_t>     m_matrixDim;
    std::vector<size_t>     m_matrixDimProducts;
        
    //The coefficient array
    //TODO: This could maybe be sped up by making it just one array all next to each other.
    std::vector<size_t>     m_arraySizes;
    std::vector<double*>    m_arrays;
    
    //For evaluation
    double                  m_currPoint;
    double                  m_currDouble;
};


#endif /* PowerBasisPolynomial_h */
