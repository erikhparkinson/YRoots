//
//  IntervalChecker2D.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/11/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalChecker2D_ipp
#define IntervalChecker2D_ipp

template <>
void IntervalChecker<Dimension::Two>::runQuadraticCheck(ChebyshevApproximation<Dimension::Two>& _approximation)
{
    enum WhichInterval {
        NegXNegY = 0,
        NegXPosY = 1,
        PosXNegY = 2,
        PosXPosY = 3
    };
    
    //If it has all been thrown out we are done.
    if(m_intervalMask[0] && m_intervalMask[1] && m_intervalMask[2] && m_intervalMask[3]) {
        return;
    }
    
    //Assume we can throw out everything remaining. Throughout the check we will try to switch
    //the bools to false. At the end if either vector is false we can throw it out.
    for(size_t i = 0; i < m_throwOutMask.size(); i++) {
        m_throwOutMask[i] = !m_intervalMask[i];
    }
    
    //The quadratic part is c3(2x^2-1)+c4xy+c5(2y^2-1)+c1x+c2y+c0.
    double* array = _approximation.getArray();
    _approximation.sumAbsValues();
    double error = _approximation.getSumAbsVal() + _approximation.getApproximationError();
    size_t degree = _approximation.getDegree();
    size_t sideLength = _approximation.getSideLength();

    //Get the variables
    double c5 = degree > 1 ? array[2*sideLength] : 0.0;
    double c4 = degree > 1 ? array[sideLength + 1] : 0.0;
    double c3 = degree > 1 ? array[2] : 0.0;
    double c2 = degree > 0 ? array[sideLength] : 0.0;
    double c1 = degree > 0 ? array[1] : 0.0;
    double c0 = array[0];

    error -= (std::abs(c5) + std::abs(c4) + std::abs(c3) + std::abs(c2) + std::abs(c1) + std::abs(c0));
    error = error < 0.0 ? 0.0 : error; //TODO: If this happens we have a problem
    
    //Repeated values from Horner's method
    double k0 = c0-c3-c5;
    double k3 = 2*c3;
    double k5 = 2*c5;
    
    //Other numbers that are repeated a bunch
    constexpr double midPoint = 2*m_randomIntervalDivider - 1;
    constexpr double midPointSqrd = midPoint*midPoint;
    constexpr double midPointChebSqrd = 2*midPointSqrd - 1;

    //Check the midpoint
    double valueMidMid = (c5+c3)*midPointChebSqrd + c4*midPointSqrd + (c2+c1)*midPoint + c0;
    EvalSign evalSignMidMid = getEvalSign(valueMidMid, error);
    
    //Check the other 4 sides
    double valueNegMid = c5*midPointChebSqrd + (c2-c4)*midPoint + c0+c3-c1;
    EvalSign evalSignNegMid = getEvalSign(valueNegMid, error);
    
    double valuePosMid = valueNegMid + 2*(c4*midPoint + c1);
    EvalSign evalSignPosMid = getEvalSign(valuePosMid, error);
    
    double valueMidNeg = c3*midPointChebSqrd + (c1-c4)*midPoint + c0+c5-c2;
    EvalSign evalSignMidNeg = getEvalSign(valueMidNeg, error);
    
    double valueMidPos = valueMidNeg + 2*(c4*midPoint + c2);
    EvalSign evalSignMidPos = getEvalSign(valueMidPos, error);

    //Check the 4 Corners
    double valuePosPos = c0 + c1 + c2 + c3 + c4 + c5;
    EvalSign evalSignPosPos = getEvalSign(valuePosPos, error);

    double valueNegPos = valuePosPos - 2*(c1+c4);
    EvalSign evalSignNegPos = getEvalSign(valueNegPos, error);

    double valuePosNeg = valuePosPos - 2*(c2+c4);
    EvalSign evalSignPosNeg = getEvalSign(valuePosNeg, error);

    double valueNegNeg = valuePosPos - 2*(c1+c2);
    EvalSign evalSignNegNeg = getEvalSign(valueNegNeg, error);
    
    m_throwOutMask[NegXNegY] = static_cast<bool>(evalSignMidMid & evalSignNegMid & evalSignMidNeg & evalSignNegNeg);
    m_throwOutMask[NegXPosY] = static_cast<bool>(evalSignMidMid & evalSignNegMid & evalSignMidPos & evalSignNegPos);
    m_throwOutMask[PosXNegY] = static_cast<bool>(evalSignMidMid & evalSignPosMid & evalSignMidNeg & evalSignPosNeg);
    m_throwOutMask[PosXPosY] = static_cast<bool>(evalSignMidMid & evalSignPosMid & evalSignMidPos & evalSignPosPos);

    if(!(m_throwOutMask[0] || m_throwOutMask[1] || m_throwOutMask[2] || m_throwOutMask[3])) {
        return;
    }
    
    //Variables for evaluating the boundaries
    double x,y;
    //Check the x constant boundaries
    if(c5 != 0) {
        double four_c5 = 4*c5;
        //When x = middleVal
        y = (-c2-c4*midPoint)/four_c5;
        if (midPoint <= y && y <= 1) {
            double value = c0 + (c1 + c4*y)*midPoint + c3*midPointChebSqrd - c5 + (k5*y+c2)*y;
            if(!static_cast<bool>(evalSignMidMid & getEvalSign(value, error))) {
                m_throwOutMask[NegXPosY] = false;
                m_throwOutMask[PosXPosY] = false;
            }
        }
        else if( -1 <= y && y <= midPoint) {
            double value = c0 + (c1 + c4*y)*midPoint + c3*midPointChebSqrd - c5 + (k5*y+c2)*y;
            if(!static_cast<bool>(evalSignMidMid & getEvalSign(value, error))) {
                m_throwOutMask[NegXNegY] = false;
                m_throwOutMask[PosXNegY] = false;
            }
        }

        //When x = -1
        y = (-c2+c4)/four_c5;
        if (midPoint <= y && y <= 1) {
            m_throwOutMask[NegXNegY] = m_throwOutMask[NegXNegY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k3 - c1 + (k5*y+c2-c4)*y, error));
        }
        else if( -1 <= y && y <= midPoint) {
            m_throwOutMask[NegXPosY] = m_throwOutMask[NegXPosY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k3 - c1 + (k5*y+c2-c4)*y, error));
        }
        
        //When x = 1
        y = (-c2-c4)/four_c5;
        if (midPoint <= y && y <= 1) {
            m_throwOutMask[PosXNegY] = m_throwOutMask[PosXNegY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k3 + c1 + (k5*y+c2+c4)*y, error));
        }
        else if(-1 <= y && y <= midPoint) {
            m_throwOutMask[PosXPosY] = m_throwOutMask[PosXPosY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k3 + c1 + (k5*y+c2+c4)*y, error));
        }
    }

    //Check the y constant boundaries
    if(c3 != 0) {
        double four_c3 = 4*c3;
        //When y = middleVal
        x = (-c1-c4*midPoint)/four_c3;
        if (midPoint <= x && x <= 1) {
            double value = c0 + (c2 + c4*x)*midPoint + c5*midPointChebSqrd - c3 + (k3*x+c1)*x;
            if(!static_cast<bool>(evalSignMidMid & getEvalSign(value, error))) {
                m_throwOutMask[PosXNegY] = false;
                m_throwOutMask[PosXPosY] = false;
            }
        }
        else if( -1 <= x && x <= midPoint) {
            double value = c0 + (c2 + c4*x)*midPoint + c5*midPointChebSqrd - c3 + (k3*x+c1)*x;
            if(!static_cast<bool>(evalSignMidMid & getEvalSign(value, error))) {
                m_throwOutMask[NegXNegY] = false;
                m_throwOutMask[NegXPosY] = false;
            }
        }
        
        //When y = -1
        x = (-c1+c4)/four_c3;
        if (midPoint <= x && x <= 1) {
            m_throwOutMask[PosXNegY] = m_throwOutMask[PosXNegY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k5 - c2 + (k3*x+c1-c4)*x, error));
        }
        else if( -1 <= x && x <= midPoint) {
            m_throwOutMask[NegXNegY] = m_throwOutMask[NegXNegY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k5 - c2 + (k3*x+c1-c4)*x, error));
        }
        
        //When y = 1
        x = (-c1-c4)/four_c3;
        if (midPoint <= x && x <= 1) {
            m_throwOutMask[PosXPosY] = m_throwOutMask[PosXPosY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k5 + c2 + (k3*x+c1+c4)*x, error));
        }
        else if( -1 <= x && x <= midPoint) {
            m_throwOutMask[NegXPosY] = m_throwOutMask[NegXPosY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + k5 + c2 + (k3*x+c1+c4)*x, error));
        }
    }

    //Check the interior
    //Comes from solving dx, dy = 0
    //Dx: 4c3x +  c4y = -c1    Matrix inverse is  [4c5  -c4]
    //Dy:  c4x + 4c5y = -c2                       [-c4  4c3]
    double det = 16 * c3 * c5 - c4*c4;
    if(det != 0) {
        x = (c2 * c4 - 4 * c1 * c5) / det;
        y = (c1 * c4 - 4 * c2 * c3) / det;

        if(midPoint <= x && x <= 1) {
            if(midPoint <= y && y <= 1) {
                m_throwOutMask[PosXPosY] = m_throwOutMask[PosXPosY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + (k3*x+c1+c4*y)*x + (k5*y+c2)*y, error));
            }
            else if(-1 <= y && y <= midPoint) {
                m_throwOutMask[PosXNegY] = m_throwOutMask[PosXNegY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + (k3*x+c1+c4*y)*x + (k5*y+c2)*y, error));
            }
        }
        else if(-1 <= x && x <= midPoint) {
            if(midPoint <= y && y <= 1) {
                m_throwOutMask[NegXPosY] = m_throwOutMask[NegXPosY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + (k3*x+c1+c4*y)*x + (k5*y+c2)*y, error));
            }
            else if(-1 <= y && y <= midPoint) {
                m_throwOutMask[NegXNegY] = m_throwOutMask[NegXNegY] && static_cast<bool>(evalSignMidMid & getEvalSign(k0 + (k3*x+c1+c4*y)*x + (k5*y+c2)*y, error));
            }
        }
    }
    
    //Update the m_intervalMask
    for(size_t i = 0; i < m_throwOutMask.size(); i++) {
        m_intervalMask[i] = m_intervalMask[i] || m_throwOutMask[i];
    }
}

template <>
double IntervalChecker<Dimension::Two>::getBoundingInterval(std::vector<ChebyshevApproximation<Dimension::Two>>& _chebyshevApproximations) {
    //Constants for this function
    /*const double MIN_MOVE = 0.001;
    const double MIN_HEIGHT = 0.0;
    const double WIDTH_MULTIPLIER = 1.01;*/

    /*
    double* array1 = _chebyshevApproximations[0].getArray();
    double* array2 = _chebyshevApproximations[1].getArray();
    size_t arraySize1 = _chebyshevApproximations[0].getPartialSideLength();
    size_t arraySize2 = _chebyshevApproximations[1].getPartialSideLength();
    size_t arrayFullSize1 = _chebyshevApproximations[0].getSideLength();
    size_t arrayFullSize2 = _chebyshevApproximations[1].getSideLength();
    double error1 = _chebyshevApproximations[0].getApproximationError();
    double error2 = _chebyshevApproximations[1].getApproximationError();

    //Updates the size of the vectors if needed
    if(unlikely(arraySize1 > m_biggestArraySizeChecked[0])) {
        m_reducedChebEvals[0].resize(arraySize1);
        m_biggestArraySizeChecked[0] = arraySize1;
    }
    if(unlikely(arraySize2 > m_biggestArraySizeChecked[1])) {
        m_reducedChebEvals[1].resize(arraySize2);
        m_biggestArraySizeChecked[1] = arraySize2;
    }
    
    //Get the Lipshitz Constants of the appoximations in both dimensions
    double lConst1X,lConst2X,lConst1Y,lConst2Y;
    getLipshitzConstant2D(array1, arraySize1, arrayFullSize1, lConst1X, lConst1Y);
    getLipshitzConstant2D(array2, arraySize2, arrayFullSize2, lConst2X, lConst2Y);

    //Set the start interval
    double a = -1;
    double b = 1;
        
    //Solve
    //Run the values
    bool makingProgress = true;
    size_t currOptDim = 0;
    while (a < b && makingProgress && currOptDim < 2) {
        makingProgress = false;
        //Get the evaluations down a dimension
        evaluateCheb2D(array1, arraySize1, arrayFullSize1, a, currOptDim, m_reducedChebEvals[0][0]);
        evaluateCheb2D(array2, arraySize2, arrayFullSize2, b, currOptDim, m_reducedChebEvals[1][0]);

        //Get the Lipshitz constants down a dimension
        double reducedLConst1 = getLipshitzConstant1D(m_reducedChebEvals[0][0].data(), arraySize1);
        double reducedLConst2 = getLipshitzConstant1D(m_reducedChebEvals[1][0].data(), arraySize2);
        
        //Determine how far we want to step and the goal heigths we need to prove
        double goalWidth = (b-a) * WIDTH_MULTIPLIER;
        double goalHeight1, goalHeight2;
        if(currOptDim == 0) {
            goalHeight1 = goalWidth*lConst1X;
            goalHeight2 = goalWidth*lConst2X;
        }
        else {
            goalHeight1 = goalWidth*lConst1Y;
            goalHeight2 = goalWidth*lConst2Y;
        }
        goalHeight1 = std::max(goalHeight1, MIN_HEIGHT);
        goalHeight2 = std::max(goalHeight2, MIN_HEIGHT);

        double a2, b2;
        a2 = 0;
        */
        
        /*

        #Try to prove the height we need
        currentVal = otherInterval[0]
        while currentVal < otherInterval[1]:
            val1 = abs(chebval1D(currentVal, M1Line)) - error1
            val2 = abs(chebval1D(currentVal, M2Line)) - error2
            while val1 < goalHeight1 and val2 < goalHeight2:
                goalWidth /= 2
                goalHeight1 /= 2
                goalHeight2 /= 2
                if goalWidth < MIN_LINE_MOVE:
                    goalWidth = 0.0
                    break
            if goalWidth == 0.0:
                break
            delta1 = (val1 - goalHeight1) / lConstLine1 if lConstLine1 != 0 else (val1 - goalHeight1)*np.inf
            delta2 = (val2 - goalHeight2) / lConstLine2 if lConstLine2 != 0 else (val2 - goalHeight2)*np.inf
            delta = max(delta1, delta2)
            if delta > 0:
                currentVal += delta
            if delta < MIN_DELTA:
                goalWidth /= 2
                goalHeight1 /= 2
                goalHeight2 /= 2
                if goalWidth < MIN_LINE_MOVE:
                    goalWidth = 0.0
                    break
        #If we weren't able to prove the bound on the segment
        if goalWidth == 0.0:
            break
        #Move the line
        if solveUpper:
            thisInterval[1] -= goalWidth
        else:
            thisInterval[0] += goalWidth
        #Stop if we have eliminated it
        if thisInterval[0] >= thisInterval[1]:
            break


        
        
    }
         */

    return std::numeric_limits<double>::max();
}


#endif /* IntervalChecker2D_ipp */
