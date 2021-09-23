//
//  Polynomial.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/25/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef Polynomial_h
#define Polynomial_h

#include <iostream>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include "Utilities/utilities.hpp"
#include "Utilities/ErrorTracker.hpp"

struct Monomial {
    std::vector<size_t> spot;
    double coeff;
    
    void clear(size_t _rank) {
        spot.clear();
        spot.resize(_rank, 0);
        coeff = 0.0;
        m_readyToEval = false;
    }
    
    template<typename ReturnType>
    ReturnType evaluate(const std::vector<double>& _values) const {
        //This is just used for the unit tests for now
        assert(spot.size() == _values.size());
        ReturnType result = coeff;
        for(size_t i = 0; i < spot.size(); i++) {
            result *= power(_values[i], spot[i]);
        }
        return result;
    }
    
    void evaluateGrid(const std::vector<std::vector<double> >& _grid, std::vector<double>& _results) {
        //Create the grid for only the dimensions we are using.
        std::vector<const std::vector<double>*> usedGrid;
        std::vector<size_t> usedSpots;
        for(size_t dim = 0; dim < _grid.size(); dim++) {
            if(m_hasDimension[dim]) {
                usedGrid.push_back(&_grid[dim]);
                usedSpots.push_back(spot[dim]);
            }
        }
        assert(usedGrid.size() != 0);
        assert(usedGrid[0]->size() != 0);
        
        //For convenience, create a reference to where the result will always be
        const size_t lastDim = m_numUsedDimensions - 1;
        double& resultValue = m_evaluationInfo[0];
        size_t resultSpot = 0;

        //Track the spot we are currently evaluating.
        std::vector<size_t> evalSpot(m_numUsedDimensions,0);

        //Evaluate the first point in the grid.
        m_evaluationInfo[lastDim] = power((*usedGrid[lastDim])[0], usedSpots[lastDim]) * coeff;
        for(size_t dim = lastDim-1; dim < m_numUsedDimensions; dim--) {
            m_evaluationInfo[dim] = power((*usedGrid[dim])[0], usedSpots[dim]) * m_evaluationInfo[dim+1];
        }
        _results[resultSpot++] = resultValue;
        
        //Iterate through all the combinations
        size_t spotToInc = 0;
        while (spotToInc < m_numUsedDimensions) {
            while(++evalSpot[spotToInc] < _grid[0].size()) {
                //Update the evaluations of the higher dimensions if needed
                while(spotToInc > 0) {
                    m_evaluationInfo[spotToInc] = power((*usedGrid[spotToInc])[evalSpot[spotToInc]], usedSpots[spotToInc]) * (spotToInc == lastDim ? coeff : m_evaluationInfo[spotToInc+1]);
                    spotToInc--;
                }
                //Evaluate the final dimension and get the result
                m_evaluationInfo[0] = power((*usedGrid[0])[evalSpot[0]], usedSpots[0]) * (0 == lastDim ? coeff : m_evaluationInfo[1]);
                _results[resultSpot++] = resultValue;
            }
            evalSpot[spotToInc] = 0;
            spotToInc++;
        }
    }
    
    void prepEvaluation() {
        m_readyToEval = true;
        //Get the dimension info
        m_hasDimension.clear();
        m_hasDimension.resize(spot.size(), false);
        m_numUsedDimensions = 0;
        for(size_t i = 0; i < m_hasDimension.size(); i++) {
            m_hasDimension[i] = spot[i] != 0;
            m_numUsedDimensions += spot[i] != 0;
        }

        //Get the evaluate grid info
        m_evaluationInfo.resize(m_numUsedDimensions, 0.0);
    }
    
    //Overwrite multiplication by a monomial
    Monomial& operator*=(const Monomial& rhs)
    {
        assert(spot.size() == rhs.spot.size());
        for(size_t i = 0; i < spot.size(); i++) {
            spot[i] += rhs.spot[i];
        }
        coeff *= rhs.coeff;
        return *this; // return the result by reference
    }
    friend Monomial operator*(Monomial lhs, const Monomial& rhs)
    {
        lhs *= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    
    const std::vector<bool>& getHasDimension() const {
        return m_hasDimension;
    }
    
private:
    bool                m_readyToEval;
    std::vector<bool>   m_hasDimension;
    size_t              m_numUsedDimensions;
    std::vector<double> m_evaluationInfo;
};

//Custom struct for sorting monomials. Sorts it [0,0],[0,1],[1,0],[1,1]
struct CustomMonomialLess {
    bool operator()(const Monomial& a, const Monomial& b) const {
        assert(a.spot.size() == b.spot.size());
        for(size_t i = 0; i < a.spot.size(); i++) {
            if(a.spot[i] < b.spot[i]) {
                return true;
            }
            else if(a.spot[i] > b.spot[i]) {
                return false;
            }
        }
        return false;
    }
} customMonomialLess;

struct PolynomialDimensionEvalInfo {
    PolynomialDimensionEvalInfo():
    m_lastP0(-1),
    m_lastP1(-1),
    m_rowStart(0)
    {}
    
    //When a vector of double, points in [stopPoints[i+1], stopPoints[i]) (stopPoints[0] will be -1) are evaluated with Horner's method
    //Then multiplied by value^powerMultiplier[i] and stored in results[i]
    //For example, x^5+3x^4-2x^3 would be stored as three points in the incoming vector, with a power multipler of 3.
    //The power multiplier allows quicker evaluations for sparse things like x^100 as opposed to looping over 100 0's.
    template<typename ReturnType, typename InputType>
    void evaluate(const std::vector<InputType>& _coeffs, const double _value) {
        for(size_t evalStep = 0; evalStep < m_results.size(); evalStep++) {
            int spot = m_breakPoints[evalStep+1];
            int end = m_breakPoints[evalStep];
            if(spot != end) {
                ReturnType result = _coeffs[spot--];
                while(spot > end) {
                    result = _coeffs[spot--] + _value * result;
                }
                if(m_powerMultiplier[evalStep] != 0) {
                    result *= power(_value, m_powerMultiplier[evalStep]);
                }
                getResults<ReturnType>()[evalStep] = result;
            }
            //else m_results[evalStep] = 0, but this will be initialized to 0 so always be that already.
        }
    }
        
    void addPoint(bool newPrev, int p0, int p1) {
        //p1 is the index in this dimension, p0, is the previous dimension, and newPrev is in anything is updated before that.
        //If newPrev is true or p0 increases, I'm on a new row
        //  For every number p0 jumps when newPrev is false, duplicate the last number in m_breakPoints, and add 0 to m_powerMultiplier and m_results
        //If newPrev is false and p0 stays the same, wait until we hit the max p1 value, then add it to breakpoints
        if(newPrev) { //This should be true on the first call
            if(m_breakPoints.size() == 0) {
                m_breakPoints.push_back(-1);
            }
            else {
                int bInc = m_lastP1 + 1 - m_rowStart;
                m_breakPoints.push_back(bInc + m_breakPoints.back());
            }
            m_powerMultiplier.push_back(p1);
            m_results.push_back(0);
            m_resultsErrorTracker.push_back(0);
            m_rowStart = p1;
        }
        else if(p0 != m_lastP0){
            int bInc = m_lastP1 + 1 - m_rowStart;
            m_breakPoints.push_back(bInc + m_breakPoints.back());
            for(int i = m_lastP0 + 1; i < p0; i++) { //The 0s we are skipping over
                m_breakPoints.push_back(m_breakPoints.back());
                m_powerMultiplier.push_back(0);
                m_results.push_back(0);
                m_resultsErrorTracker.push_back(0);
            }
            m_powerMultiplier.push_back(p1);
            m_results.push_back(0);
            m_resultsErrorTracker.push_back(0);
            m_rowStart = p1;
        }
        m_lastP0 = p0;
        m_lastP1 = p1;
    }
    
    void endInit() {
        int bInc = m_lastP1 + 1 - m_rowStart;
        m_breakPoints.push_back(bInc + m_breakPoints.back());
    }
    
    template<typename ReturnType>
    inline std::vector<ReturnType>& getResults() {
        
    }

    template<>
    inline std::vector<double>& getResults() {
        return m_results;
    }

    template<>
    inline std::vector<ErrorTracker>& getResults() {
        return m_resultsErrorTracker;
    }

    //m_breakPoints starts with -1 and is exactly 1 bigger than the other two.
    std::vector<int>        m_breakPoints;
    std::vector<size_t>     m_powerMultiplier;
    
    std::vector<double>         m_results;
    std::vector<ErrorTracker>   m_resultsErrorTracker;

    //For addPoints
    int                     m_lastP0;
    int                     m_lastP1;
    int                     m_rowStart;
};


class Polynomial {
public:
    Polynomial() :
    m_readyToEval(false)
    {}
    
    void setRank(size_t _rank) {
        m_rank = _rank;
    }

    double evaluateSlow(const std::vector<double>& _values) { //Only use for testing and verifying
        double result = 0;
        for(const Monomial& m : m_monomials) {
            result += m.evaluate<double>(_values);
        }
        return result;
    }

    template<typename ReturnType>
    ReturnType evaluate(const std::vector<double>& _values) {
        if(unlikely(!m_readyToEval)) {
            prepEvaluation();
        }
        if(unlikely(m_numUsedDimensions == 0)) { //Check the case of an empty polynomial.
            return static_cast<ReturnType>(m_constantTerm);
        }
        //Get the dimension to run and match it to the right value dimension
        const size_t lastDim = m_evaluationInfo.size() - 1;
        size_t valueEvalSpot = _values.size() - 1;
        while(!m_hasDimension[valueEvalSpot]) {
            valueEvalSpot--;
        }
        //Evaluate the first dimension
        m_evaluationInfo[lastDim].evaluate<ReturnType>(m_coeffs, _values[valueEvalSpot]);
        //Loop through the rest of the dimensions.
        for(size_t dim = lastDim-1; dim < m_evaluationInfo.size(); dim--) {
            while(!m_hasDimension[--valueEvalSpot]) {} //Get to the next value dimension
            m_evaluationInfo[dim].evaluate<ReturnType>(m_evaluationInfo[dim+1].getResults<ReturnType>(), _values[valueEvalSpot]);
        }
        return m_evaluationInfo[0].getResults<ReturnType>()[0];
    }
    
    void evaluateGrid(const std::vector<std::vector<double> >& _grid, std::vector<double>& _results) {
        if(unlikely(!m_readyToEval)) {
            prepEvaluation();
        }
        if(unlikely(m_numUsedDimensions == 0)) { //Check the case of an empty polynomial.
            std::fill(_results.begin(), _results.end(), m_constantTerm);
            return;
        }
        //Create the grid for only the dimensions we are using.
        std::vector<const std::vector<double>*> usedGrid;
        for(size_t dim = 0; dim < _grid.size(); dim++) {
            if(m_hasDimension[dim]) {
                usedGrid.push_back(&_grid[dim]);
            }
        }
        assert(usedGrid.size() != 0);
        assert(usedGrid[0]->size() != 0);
        
        //For convenience, create a reference to where the result will always be
        const size_t lastDim = m_numUsedDimensions - 1;
        double& resultValue = m_evaluationInfo[0].m_results[0];
        size_t resultSpot = 0;

        //Track the spot we are currently evaluating.
        std::vector<size_t> evalSpot(m_numUsedDimensions,0);

        //Evaluate the first point in the grid.
        m_evaluationInfo[lastDim].evaluate<double>(m_coeffs, (*usedGrid[lastDim])[0]);
        for(size_t dim = lastDim-1; dim < m_numUsedDimensions; dim--) {
            m_evaluationInfo[dim].evaluate<double>(m_evaluationInfo[dim+1].m_results, (*usedGrid[dim])[0]);
        }
        _results[resultSpot++] = resultValue;
        
        //Iterate through all the combinations
        size_t spotToInc = 0;
        while (spotToInc < m_numUsedDimensions) {
            while(++evalSpot[spotToInc] < _grid[0].size()) {
                //Update the evaluations of the higher dimensions if needed
                while(spotToInc > 0) {
                    m_evaluationInfo[spotToInc].evaluate<double>(spotToInc == lastDim ? m_coeffs : m_evaluationInfo[spotToInc+1].m_results, (*usedGrid[spotToInc])[evalSpot[spotToInc]]);
                    spotToInc--;
                }
                //Evaluate the final dimension and get the result
                m_evaluationInfo[0].evaluate<double>(0 == lastDim ? m_coeffs : m_evaluationInfo[1].m_results, (*usedGrid[0])[evalSpot[0]]);
                _results[resultSpot++] = resultValue;
            }
            evalSpot[spotToInc] = 0;
            spotToInc++;
        }
    }
    
    void prepEvaluation() {
        std::vector<Monomial> monomialVector;
        for(const Monomial& m : m_monomials) {
            monomialVector.push_back(m);
        }
        prepEvaluation(monomialVector);
    }
        
    void multiplyMonomial(const Monomial& _newMonomial) {
        //Pulls everything of the set and into a vector, and the pushes the multiples back in.
        std::vector<Monomial> monomialVector;
        for(const Monomial& m : m_monomials) {
            monomialVector.push_back(m);
        }
        m_monomials.clear();
        for(const Monomial& m : monomialVector) {
            m_monomials.insert(m * _newMonomial);
        }
    }
    
    void addMonomials(const std::vector<Monomial>& _newMonomials) {
        for(const Monomial& m : _newMonomials) {
            addMonomial(m);
        }
    }

    void addMonomials(const std::set<Monomial, CustomMonomialLess>& _newMonomials) {
        for(const Monomial& m : _newMonomials) {
            addMonomial(m);
        }
    }

    void addMonomial(const Monomial& _newMonomial) {
        std::set<Monomial, CustomMonomialLess>::iterator found = m_monomials.find(_newMonomial);
        if(found != m_monomials.end()) {
            Monomial toInsert = _newMonomial;
            toInsert.coeff += found->coeff;
            m_monomials.erase(found);
            m_monomials.insert(toInsert);
        }
        else {
            m_monomials.insert(_newMonomial);
        }
    }
    
    const std::set<Monomial, CustomMonomialLess>& getMonomials() const {
        return m_monomials;
    }
    
    const std::vector<bool>& getHasDimension() const {
        return m_hasDimension;
    }

    size_t getNumUsedDimensions() const {
        return m_numUsedDimensions;
    }
    
    void clear() {
        m_coeffs.clear();
        m_evaluationInfo.clear();
        m_hasDimension.clear();
        m_numUsedDimensions = 0;
        m_constantTerm = 0;
        m_monomials.clear();
    }
    
private:
    void prepEvaluation(std::vector<Monomial>& _monomials) { //Monomials gets destroyed here
        m_readyToEval = true;
        if(unlikely(_monomials.size() == 0)) { //Check the case of an empty polynomial.
            m_numUsedDimensions = 0;
            m_constantTerm = 0.0;
            return;
        }
        assert(_monomials[0].spot.size() > 0);
        
        //Sort the monomials.
        std::sort(_monomials.begin(), _monomials.end(), customMonomialLess);
        
        //Find out what dimensions are used.
        m_hasDimension.resize(m_rank, false);
        for(size_t i = 0; i < _monomials.size(); i++) {
            for(size_t j = 0; j < m_rank; j++) {
                m_hasDimension[j] = m_hasDimension[j] || _monomials[i].spot[j] > 0;
            }
        }
        m_numUsedDimensions = 0;
        for(size_t i = 0; i < m_rank; i++) {
            m_numUsedDimensions += m_hasDimension[i];
        }
        
        //Check the case of just the constant term
        if(unlikely(_monomials.size() == 1 && m_numUsedDimensions == 0)) { //Check the case of an empty polynomial.
            m_constantTerm = _monomials[0].coeff;
            return;
        }
        
        //Remove anything unused from the monomials
        for(size_t i = 0; i < _monomials.size(); i++) {
            for(size_t j = m_rank-1; j < m_rank; j--) {
                if(!m_hasDimension[j]) {
                    _monomials[i].spot.erase(_monomials[i].spot.begin()+j);
                }
            }
        }
        
        //Create the m_evaluationInfos
        m_evaluationInfo.clear();
        m_coeffs.clear();
        m_evaluationInfo.resize(m_numUsedDimensions);
        
        size_t lastDimension = m_numUsedDimensions - 1;
        size_t lastSpot = 0;
        bool prev1Change = true; //Track if a monomial of dimension less then the current one has changed.
        bool prev2Change = true; //Track if a monomial of dimension at least 2 less then the current one has changed.
        for(size_t i = 0; i < _monomials.size(); i++) {
            for(size_t j = 0; j < m_numUsedDimensions; j++) {
                //Create the evaluationInfo
                if(j == 0) {
                    m_evaluationInfo[j].addPoint(prev2Change, 0, _monomials[i].spot[j]);
                }
                else {
                    m_evaluationInfo[j].addPoint(prev2Change, _monomials[i].spot[j-1], _monomials[i].spot[j]);
                }
                
                //Set up the coeff vector
                if(j == lastDimension) {
                    if(prev1Change) {
                        m_coeffs.push_back(_monomials[i].coeff);
                    }
                    else {
                        for(size_t k = lastSpot + 1; k < _monomials[i].spot[j]; k++) {
                            m_coeffs.push_back(0);
                        }
                        m_coeffs.push_back(_monomials[i].coeff);
                    }
                    lastSpot = _monomials[i].spot[j];
                }
                
                //Record if the dimensions more than 1 before have changed
                if(i > 0) {
                    if(j > 0) {
                        prev2Change |= _monomials[i].spot[j-1] != _monomials[i-1].spot[j-1];
                    }
                    prev1Change |= _monomials[i].spot[j] != _monomials[i-1].spot[j];
                }
            }
            prev2Change = false;
            prev1Change = false;
        }
        
        for(size_t i = 0; i < m_numUsedDimensions; i++) {
            m_evaluationInfo[i].endInit();
        }
        
        //TODO: Think about: Another option to maybe simplify things is if we have a bunch of zeros between things (x^8+x^2)
        //Store how many zeros are in between each loop, and then have a switch statement that multiplies by the power, just breaks if it's zero?
        //This might be good for not very dense things, but if it's dense it's a waste having the switch every time.
        //Maybe I should have multiple eval functions, depending on if it's dense or not?
        //  Having that on every row would mean a switch statement on each of those which wouldn't be too slow, but all that could
        //  hurt the branch prediction and have other issues.
        //      I guess each dimension would be evaluated differently, maybe it's sparse in the first dimension, but then dense in the final dimension?
        
        //TODO: Check what the most efficient evaluation method would be, between what it currently has, maybe the above method, and then
        //just summing up the monomials for really sparse stuff. Prep the eval for all of them, time it quick, and then have a switch statement
        //to determine which was to actually evaluate it?
        //  It's possible it would depend on the size of the grid, but maybe not.
    }
    
private:
    //All the coefficients
    std::vector<double>                         m_coeffs;
    //The information for reducing the evaluation down each dimension.
    std::vector<PolynomialDimensionEvalInfo>    m_evaluationInfo;
    //Which dimensions are used in the polynomials
    std::vector<bool>                           m_hasDimension;
    size_t                                      m_numUsedDimensions;
    //The rank of the space we are in, size of the m_hasDimension vector.
    size_t                                      m_rank;
    //For the special case of just a constant term.
    double                                      m_constantTerm;
    
    //The monomials the make up the polynomial
    std::set<Monomial, CustomMonomialLess>      m_monomials;
    
    bool                                        m_readyToEval;
};

#endif /* Polynomial_h */
