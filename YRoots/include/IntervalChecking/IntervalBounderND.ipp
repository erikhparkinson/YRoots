//
//  IntervalBounderND.ipp
//  YRoots
//
//  Created by Erik Hales Parkinson on 8/14/21.
//  Copyright © 2021 Erik Hales Parkinson. All rights reserved.
//

#ifndef IntervalBounderND_h
#define IntervalBounderND_h

template <int Rank>
IntervalBounder<Rank>::IntervalBounder(size_t _rank):
m_rank(_rank),
m_linearTerms(m_rank, m_rank),
m_constantTerms(m_rank),
m_errorTerms(m_rank),
m_errorTermsOfLinear(m_rank),
m_rightHandSideOfLinearSystemErrors(m_rank, power(2, m_rank))
{
    m_boundingInterval.lowerBounds.resize(m_rank, 0.0);
    m_boundingInterval.upperBounds.resize(m_rank, 0.0);
    
    m_timer.registerTimer(m_timerLinearErrorSolveIndex, "Linear Error Solve");
    m_timer.registerTimer(m_timerPreconditionPolysIndex, "Precondition Polys");
    m_timer.registerTimer(m_timerLipshitzSolveIndex, "Lipshitz Solve");
    m_timer.registerTimer(m_timerChebValReduceIndex, "Cheb Val Reduce");
}

template<int Rank>
double IntervalBounder<Rank>::updateBoundingIntervalLinearErrorSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations, const std::vector<bool>& _allowedToReduceDim) {
    //Set up m_linearTerms, m_constantTerms and m_errorTerms
    for(size_t i = 0; i < m_rank; i++) { //Iterate through each approximation
        size_t spot = 1;
        size_t sideLength = _chebyshevApproximations[i].getSideLength();
        double linearTermAbsSum = 0.0;
        for(size_t j = 0; j < m_rank; j++) { //Iterate through the linear terms in the approximation
            m_linearTerms(i,j) = _chebyshevApproximations[i].getArray()[spot];
            linearTermAbsSum += std::abs(_chebyshevApproximations[i].getArray()[spot]);
            spot *= sideLength;
        }
        //Set the constant term
        m_constantTerms(i) = _chebyshevApproximations[i].getArray()[0];
        //Set the error terms
        m_errorTerms(i) = _chebyshevApproximations[i].getApproximationError();
        //Add everything but the linear stuff to the error terms
        m_errorTermsOfLinear(i) = (_chebyshevApproximations[i].getSumAbsVal() - linearTermAbsSum - std::abs(m_constantTerms(i))) + m_errorTerms(i);
    }
    
    //Get the matrix with the errors. TODO: Make this faster/do ahead of time computations
    size_t power2 = power(2, m_rank);
    for(size_t i = 0; i < m_rank; i++) {
        for(size_t j = 0; j < power2; j++) {
            if((j >> i) % 2 == 1) {
                m_rightHandSideOfLinearSystemErrors(i,j) = -m_constantTerms(i) + m_errorTermsOfLinear(i);
            }
            else {
                m_rightHandSideOfLinearSystemErrors(i,j) = -m_constantTerms(i) - m_errorTermsOfLinear(i);
            }
        }
    }
    
    //Solve the linear terms with errors.
    m_linearTermsQR = m_linearTerms.colPivHouseholderQr(); //TODO: Check if this is singular!!! Scale everything by the inf norms?
    m_linearTermsQR.setThreshold(1e-10);
    if(!m_linearTermsQR.isInvertible()) {
        return m_boundingInterval.getArea();
    }
    m_linearSystemWithErrorResult = m_linearTermsQR.solve(m_rightHandSideOfLinearSystemErrors);
    
    //Get the bounds from the linear system. //TODO: Use Eigen's fast minCoeff Call?
    for(size_t i = 0; i < m_rank; i++) {
        if(!_allowedToReduceDim[i] || true) { //TODO: Remove after testing
            continue;
        }
        double min = 1;
        double max = -1;
        for(size_t j = 0; j < power2; j++) {
            min = std::min(min, m_linearSystemWithErrorResult(i,j));
            max = std::max(max, m_linearSystemWithErrorResult(i,j));
        }
        m_boundingInterval.lowerBounds[i] = std::max(m_boundingInterval.lowerBounds[i], min);
        m_boundingInterval.upperBounds[i] = std::min(m_boundingInterval.upperBounds[i], max);
        if(m_boundingInterval.lowerBounds[i] >= m_boundingInterval.upperBounds[i]) {
            return 0.0; //Throw out the interval!
        }
    }
    
    //Return the resulting area
    return m_boundingInterval.getArea();
}

template <int Rank>
void IntervalBounder<Rank>::preconditionPolynomials(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations) { //TODO: Speed this up!
    //Create the polynomials matrix.
    size_t maxPartialSideLength = 0;
    for (size_t i = 0 ; i < _chebyshevApproximations.size(); i++) {
        maxPartialSideLength = std::max(maxPartialSideLength, _chebyshevApproximations[i].getPartialSideLength());
    }
    Eigen::MatrixXd polynomials(m_rank, power(maxPartialSideLength, m_rank));
    m_preconditionPolysDegree = maxPartialSideLength;
    
    //Initialize the polynomials matrix.
    for(size_t polyNum = 0; polyNum < _chebyshevApproximations.size(); polyNum++) {
        const size_t currSideLength = _chebyshevApproximations[polyNum].getSideLength();
        const size_t currPartialSideLength = _chebyshevApproximations[polyNum].getPartialSideLength();
        const double* approx = _chebyshevApproximations[polyNum].getArray();
        
        //Set up the needed variables
        std::vector<size_t> inputSpot(m_rank,0);
        std::vector<size_t> multipliers1(m_rank, 1); //Spots in the new vector
        std::vector<size_t> multipliers2(m_rank, 1); //Spots in thd old vector
        for(size_t i = 1; i < m_rank; i++) {
            multipliers1[i] = maxPartialSideLength*multipliers1[i-1];
            multipliers2[i] = currSideLength*multipliers2[i-1];
        }
        
        //Iterate through all the combinations
        size_t spotToInc = 0;
        polynomials(polyNum, 0) = approx[0];
        while (spotToInc < m_rank) {
            while(++inputSpot[spotToInc] < maxPartialSideLength) {
                size_t spot1 = 0;
                size_t spot2 = 0;
                bool validSpot = true;
                for (size_t i = 0; i < m_rank; i++) {
                    spot1 += inputSpot[i]*multipliers1[i];
                    spot2 += inputSpot[i]*multipliers2[i];
                    validSpot &= inputSpot[i] < currPartialSideLength;
                }
                polynomials(polyNum, spot1) = validSpot ? approx[spot2] : 0.0;
                
                if(spotToInc != 0) {
                    spotToInc = 0;
                }
            }
            inputSpot[spotToInc] = 0;
            spotToInc++;
        }
    }
    
    //Get the preconditioned polynomials
    if(m_linearTermsQR.isInvertible()) {
        m_linearTermsInverse = m_linearTerms.inverse(); //TODO: Figure out how to not invert the matrix here!
        m_preconditionedErrors = m_linearTermsInverse.cwiseAbs() * m_errorTerms;
        m_preconditionedPolys = m_linearTermsQR.solve(polynomials);
    }
    else {
        m_preconditionPolysDegree = 0;
    }
    
    /*std::cout<<"Polys:\n";
    printMatrix(polynomials);
    std::cout<<"Linear Terms:\n";
    printMatrix(m_linearTerms);
    std::cout<<"Preconditioned Polys:\n";
    printMatrix(m_preconditionedPolys);
    std::cout<<"Errors:\n";
    printMatrix(m_preconditionedErrors);*/
}

template <int Rank>
double IntervalBounder<Rank>::computeLipshitzConstant(const Eigen::VectorXd& poly, size_t dim) {
    //poly is of length m_preconditionPolysDegree^m_rank
    //Sum the squares of the spot times the place mod (m_preconditionPolysDegree)^(dim+1) or something like that.
    double value = 0.0;
    const size_t mod = power(m_preconditionPolysDegree, dim);
    for(int spot = 0; spot < poly.size(); spot++) {
        const size_t n = (spot / mod) % m_preconditionPolysDegree;
        value += std::abs(poly(spot)) * n * n;
    }
    return value;
}

template <int Rank>
void IntervalBounder<Rank>::chebValReduce(const Eigen::VectorXd& poly, Eigen::VectorXd& result, size_t dim, double value) {
    //Evaluates poly in dimension dim at value and returns the polynomials reduced down a dimension.
    size_t jumpBetweenPoints = power(m_preconditionPolysDegree, dim);
    size_t secondJump = power(m_preconditionPolysDegree, m_rank - dim - 1);
    
    //Loop through sets of m_preconditionPolysDegree points that are jumpBetweenPoints apart.
    //Start point starts at 0, then increases by 1 for a count of jumpBetweenPoints. Then increases by jumpBetweenPoints*m_preconditionPolysDegree
    
    size_t spot = 0;
    size_t resultSpot = 0;
    const double valueTimes2 = 2*value;
    for(size_t i = 0; i < secondJump; i++) {
        for(size_t j = 0; j < jumpBetweenPoints; j++) {
            //Evaluate the poly at point spot + m_preconditionPolysDegree jumpBetweenPoints. Store in resultSpot
            if(m_preconditionPolysDegree == 1) {
                result(resultSpot) = poly(spot);
            }
            else {
                double c1, c0, temp;
                size_t index = spot + (m_preconditionPolysDegree - 1) * jumpBetweenPoints;
                c1 = poly(index);
                index -= jumpBetweenPoints;
                c0 = poly(index);
                //Clenshaw recursion
                for(size_t k = 2; k < m_preconditionPolysDegree; k++) {
                    index -= jumpBetweenPoints;
                    temp = c0;
                    c0 = poly(index) - c1;
                    c1 = temp + c1*valueTimes2;
                }
                result(resultSpot) = c0 + c1*value;
            }
            
            //Eval spot
            spot++;
            resultSpot++;
        }
        spot += jumpBetweenPoints*(m_preconditionPolysDegree-1);
    }
}

template<>
double IntervalBounder<2>::getExtremeAbsVal(const Eigen::VectorXd& poly, const Interval& boundingInterval, size_t dim) {
    //Do a cubic check for the 2D case
    if(poly.size() > 3) {
        return optimizeLine3(poly, boundingInterval.lowerBounds[1-dim], boundingInterval.upperBounds[1-dim]);
    }
    else if (poly.size() > 2) {
        return optimizeLine2(poly, boundingInterval.lowerBounds[1-dim], boundingInterval.upperBounds[1-dim]);
    }
    else if (poly.size() > 1) {
        return optimizeLine1(poly, boundingInterval.lowerBounds[1-dim], boundingInterval.upperBounds[1-dim]);
    }
    else {
        return optimizeLine0(poly);
    }
}

template<int Rank>
double IntervalBounder<Rank>::getExtremeAbsVal(const Eigen::VectorXd& poly, const Interval& boundingInterval, size_t dim) {
    //TODO: Think about making this an ND Quadratic check.
    double value = std::abs(poly(0));
    //Loop over everything
    for(int i = 1; i < poly.size(); i++) {
        //Bound anything with linear terms by the max's of the linear.
        double multiplier = 1.0;
        size_t tempVal = i;
        for(size_t d = 0; d < m_rank; d++) {
            if(d == dim) {
                continue;
            }
            if(tempVal % m_preconditionPolysDegree == 1) {
                multiplier *= std::max(std::abs(boundingInterval.lowerBounds[d]), std::abs(boundingInterval.upperBounds[d]));
            }
            tempVal /= m_preconditionPolysDegree;
        }
        //Subtract the coefficient
        value -= multiplier * std::abs(poly(i));
    }
        
    return value > 0.0 ? value : 0.0;
}

template<int Rank>
double IntervalBounder<Rank>::getLiphsitzBoundIncreaseND(const Eigen::VectorXd& poly, double error, double lipshitzConstant, const Interval& boundingInterval, size_t dim) {
    double extremeAbsVal = getExtremeAbsVal(poly, boundingInterval, dim);
    if(extremeAbsVal > error) {
        return (extremeAbsVal - error) / lipshitzConstant;
    }
    else {
        return 0.0;
    }
}


template <int Rank>
double IntervalBounder<Rank>::updateBoundingIntervalLipshitzSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations, const std::vector<bool>& _allowedToReduceDim) {
    if(m_preconditionPolysDegree == 0) {
        return m_boundingInterval.getArea();
    }
    //Get the Lipshitz Constants of the functions.
    Eigen::VectorXd lipshitzConstants(m_rank);
    for(size_t i = 0; i < m_rank; i++) {
        lipshitzConstants(i) = computeLipshitzConstant(m_preconditionedPolys.row(i), i);
    }
    
    //Create the reduced polynomials
    const size_t reducedPolySize = power(m_preconditionPolysDegree, m_rank-1);
    Eigen::VectorXd reducedPoly1(reducedPolySize);
    Eigen::VectorXd reducedPoly2(reducedPolySize);
    
    static const double MIN_MOVE = 1e-3;
    bool changed = true;
    while(changed) {
        changed = false;
        for(size_t dim = 0; dim < m_rank; dim++) {
            if(!_allowedToReduceDim[dim]) {
                continue;
            }
            //Evaluate the poly at the endpoints to get it down a dimension.
            m_timer.startTimer(m_timerChebValReduceIndex);
            chebValReduce(m_preconditionedPolys.row(dim), reducedPoly1, dim, m_boundingInterval.lowerBounds[dim]);
            chebValReduce(m_preconditionedPolys.row(dim), reducedPoly2, dim, m_boundingInterval.upperBounds[dim]);
            m_timer.stopTimer(m_timerChebValReduceIndex);

            //Call get Bound Increase on each function
            double delta1 = getLiphsitzBoundIncreaseND(reducedPoly1, m_preconditionedErrors(dim), lipshitzConstants(dim), m_boundingInterval, dim);
            changed |= delta1 > MIN_MOVE;
            m_boundingInterval.lowerBounds[dim] += delta1;

            double delta2 = getLiphsitzBoundIncreaseND(reducedPoly2, m_preconditionedErrors(dim), lipshitzConstants(dim), m_boundingInterval, dim);
            changed |= delta2 > MIN_MOVE;
            m_boundingInterval.upperBounds[dim] -= delta2;
            
            if(m_boundingInterval.lowerBounds[dim] >= m_boundingInterval.upperBounds[dim]) {
                return 0.0;
            }
        }
    }
    return m_boundingInterval.getArea();
}

template <int Rank>
double IntervalBounder<Rank>::computeBoundingInterval(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations, const std::vector<bool>& _allowedToReduceDim) {
    //Reset the bounding interval to [-1,1]
    for(size_t i = 0; i < m_rank; i++) {
        m_boundingInterval.lowerBounds[i] = -1.0;
        m_boundingInterval.upperBounds[i] = 1.0;
    }
    m_boundingInterval.areaFound = false;
    //Try the linear error solve
    m_timer.startTimer(m_timerLinearErrorSolveIndex);
    double size = updateBoundingIntervalLinearErrorSolve(_chebyshevApproximations, _allowedToReduceDim);
    m_timer.stopTimer(m_timerLinearErrorSolveIndex);

    if(size == 0.0) {
        return 0.0;
    }
    //return m_boundingInterval.getArea(); //TODO: Remove this once the line check is in.
    m_boundingInterval.areaFound = false;
    
    //TODO: Don't call the constant term check, as the updateBoundingIntervalLinearErrorSolve is better
    //Figure out how to do the absValueSum better.

    
    //Call the Preconditioner.
    m_timer.startTimer(m_timerPreconditionPolysIndex);
    preconditionPolynomials(_chebyshevApproximations);
    m_timer.stopTimer(m_timerPreconditionPolysIndex);

    //Call the Lipshitz Solve
    
    //TODO: Should I just return without the harder checks if it's already be narrowed down to a small enough interval?
    
    //Run the Lipshitz Solve
    m_timer.startTimer(m_timerLipshitzSolveIndex);
    size = updateBoundingIntervalLipshitzSolve(_chebyshevApproximations, _allowedToReduceDim);
    m_timer.stopTimer(m_timerLipshitzSolveIndex);

    return size;
}

#endif /* IntervalBounderND_h */
