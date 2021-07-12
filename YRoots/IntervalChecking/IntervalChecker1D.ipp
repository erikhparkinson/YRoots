//
//  IntervalChecker.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/7/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalChecker1D_ipp
#define IntervalChecker1D_ipp

template <>
void IntervalChecker<Dimension::One>::runQuadraticCheck(ChebyshevApproximation<Dimension::One>& _approximation)
{
    //The quadratic part is a(2x^2-1)+bx+c
    double* array = _approximation.getArray();
    _approximation.sumAbsValues();
    double error = _approximation.getSumAbsVal() + _approximation.getApproximationError();
    size_t degree = _approximation.getDegree();
    
    //Get the variables
    double a = degree > 1 ? array[2] : 0.0;
    double b = degree > 0 ? array[1] : 0.0;
    double c = array[0];
    error -= (std::abs(a) + std::abs(b) + std::abs(c));
    error = error < 0.0 ? 0.0 : error; //TODO: If this happens we have a problem
    
    //Check the midpoint
    double midPoint = 2*m_randomIntervalDivider - 1;
    double valueMid = a * (2*midPoint*midPoint - 1) + b * midPoint + c;
    double valuePos = a + b + c;
    double valueNeg = a - b + c;
    
    EvalSign evalSignMid = getEvalSign(valueMid, error);
    EvalSign evalSignPos = getEvalSign(valuePos, error);
    EvalSign evalSignNeg = getEvalSign(valueNeg, error);

    m_intervalMask[0] = (evalSignMid & evalSignNeg);
    m_intervalMask[1] = (evalSignMid & evalSignPos);
    
    //Take the min into account
    if(a != 0) {
        //Min occurs when 4ax + b = 0
        double minValSpot = -b/(4*a);
        //Evaluating gives a(2*b^2/16a^2 - 1) + b(-b/4a) + c = b^2/8a - b^2/4a + c - a = - b^2/8a + c - a
        double minVal = -b*b/(8*a) + c - a;
        int minValSign = getEvalSign(minVal, error);
        bool newBool = (evalSignMid & minValSign);
        if(minValSpot > -1 && minValSpot < 1) {
            if(minValSpot < midPoint) {
                m_intervalMask[0] = m_intervalMask[0] && newBool;
            }
            else {
                m_intervalMask[1] = m_intervalMask[1] && newBool;
            }
        }
    }
}

template <>
double IntervalChecker<Dimension::One>::getBoundingInterval(std::vector<ChebyshevApproximation<Dimension::One> >& _chebyshevApproximations) {
    return std::numeric_limits<double>::max();
    
    //Constants for this function
    const double MIN_MOVE = 0.001;

    double* array = _chebyshevApproximations[0].getArray();
    size_t arraySize = _chebyshevApproximations[0].getPartialSideLength();
    
    //Get the Lipshitz constant
    double lipshitzConstant = getLipshitzConstant1D(array, arraySize);
    
    double error = _chebyshevApproximations[0].getApproximationError();
    double a = -1;
    double b = 1;
    
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

    //Return the results
    if(a > b) {
        return 0.0;
    }
    else {
        m_boundingInterval.lowerBounds[0] = a;
        m_boundingInterval.upperBounds[0] = b;
        return b - a;
    }
}

#endif /* IntervalChecker1D_ipp */
