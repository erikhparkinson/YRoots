//
//  IntervalBounder1D.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 8/14/21.
//  Copyright © 2021 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalBounder1D_h
#define IntervalBounder1D_h

template<>
double IntervalBounder<1>::updateBoundingIntervalLinearErrorSolve(std::vector<ChebyshevApproximation<1> >& _chebyshevApproximations) {
    //The Constant array information
    const double* array = _chebyshevApproximations[0].getArray();
    const size_t arraySize = _chebyshevApproximations[0].getPartialSideLength();
    const double error = _chebyshevApproximations[0].getApproximationError();

    //The the error of the higher than linear terms
    double postLinearError = 0.0;
    for (size_t i = 2; i < arraySize; i++) {
        postLinearError += std::abs(array[i]);
    }
    
    //A Linear system ax + b +- e has solutions in [-b/a +- e/a]
    double linearSystemCenter = array[0] / array[1];
    double linearSystemWidth = std::abs((error + postLinearError)/ array[1]);

    //Update the Boudning Interval with the bounds from solving the linear system.
    m_boundingInterval.lowerBounds[0] = std::max(m_boundingInterval.lowerBounds[0], linearSystemCenter - linearSystemWidth);
    m_boundingInterval.upperBounds[0] = std::max(m_boundingInterval.upperBounds[0], linearSystemCenter + linearSystemWidth);
    
    return m_boundingInterval.upperBounds[0] - m_boundingInterval.lowerBounds[0];
}

template <>
double IntervalBounder<1>::updateBoundingIntervalLipshitzSolve(std::vector<ChebyshevApproximation<1> >& _chebyshevApproximations) {
    //Constants for this function
    const double MIN_MOVE = 1e-3;

    //The Constant array information
    const double* array = _chebyshevApproximations[0].getArray();
    const size_t arraySize = _chebyshevApproximations[0].getPartialSideLength();
    const double error = _chebyshevApproximations[0].getApproximationError();

    //Get the Lipshitz constant
    double lipshitzConstant = getLipshitzConstant1D(array, arraySize);
        
    //Write the bounds as a,b for simplicity
    double& a = m_boundingInterval.lowerBounds[0];
    double& b = m_boundingInterval.upperBounds[0];

    //Get the intitial evaluations
    double evalA = std::abs(evaluateCheb1D(array, arraySize, a));
    double evalB = std::abs(evaluateCheb1D(array, arraySize, b));

    //If the lipshitzConstant is 0, then we can throw it out if we get an eval above the error.
    if(lipshitzConstant == 0) {
        if(evalA > error || evalB > error) {
            return 0.0;
        }
        return std::numeric_limits<double>::max();
    }
    
    //Run the values
    bool makingProgress = true;
    while (a < b && makingProgress) {
        makingProgress = false;
        
        if(evalA > error) {
            double move = (evalA - error) / lipshitzConstant;
            a += move;
            makingProgress |= (move > MIN_MOVE);
            evalA = std::abs(evaluateCheb1D(array, arraySize, a));
        }
        if(evalB > error) {
            double move = (evalB - error) / lipshitzConstant;
            b -= move;
            makingProgress |= (move > MIN_MOVE);
            evalB = std::abs(evaluateCheb1D(array, arraySize, b));
        }
    }

    //Return the area remaining.
    if(a > b) {
        return 0.0;
    }
    else {
        return m_boundingInterval.upperBounds[0] - m_boundingInterval.lowerBounds[0];
    }
}

template <>
double IntervalBounder<1>::computeBoundingInterval(std::vector<ChebyshevApproximation<1> >& _chebyshevApproximations) {
    m_boundingInterval.areaFound = false;
    if(updateBoundingIntervalLinearErrorSolve(_chebyshevApproximations) == 0) {
        return 0.0;
    }
    m_boundingInterval.areaFound = false;
    return updateBoundingIntervalLipshitzSolve(_chebyshevApproximations);
}

#endif /* IntervalBounder1D_h */
