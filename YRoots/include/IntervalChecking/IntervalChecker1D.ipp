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
void IntervalChecker<1>::runQuadraticCheck(ChebyshevApproximation<1>& _approximation)
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

#endif /* IntervalChecker1D_ipp */
