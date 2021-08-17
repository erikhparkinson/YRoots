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
m_rightHandSideOfLinearSystemErrors(m_rank, 2)
{
    m_boundingInterval.lowerBounds.resize(m_rank, 0.0);
    m_boundingInterval.upperBounds.resize(m_rank, 0.0);
}

template<int Rank>
double IntervalBounder<Rank>::updateBoundingIntervalLinearErrorSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations) {
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
        m_constantTerms(i) = -_chebyshevApproximations[i].getArray()[0];
        //Set the error terms
        m_errorTerms(i) = _chebyshevApproximations[i].getApproximationError();
        //Add everything but the linear stuff to the error terms
        m_errorTermsOfLinear(i) = (_chebyshevApproximations[i].getSumAbsVal() - linearTermAbsSum - m_constantTerms(i)) + m_errorTerms(i);
    }
    
    //Get the matrix with the errors. TODO: Precompute some of these?
    size_t powerLimit = power(2, m_rank);
    for(size_t i = 0; i < m_rank; i++) {
        m_rightHandSideOfLinearSystemErrors(i,0) = -m_constantTerms(i);
        m_rightHandSideOfLinearSystemErrors(i,1) = m_errorTermsOfLinear(i);
    }
    
    //Solve the linear terms with errors.
    m_linearTermsQR = m_linearTerms.colPivHouseholderQr(); //TODO: Check if this is singular!!! Scale everything by the inf norms?
    m_linearSystemWithErrorResult = m_linearTermsQR.solve(m_rightHandSideOfLinearSystemErrors);
    
    //Get the bounds from the linear system. //TODO: Use Eigen's fast minCoeff Call?
    for(size_t i = 0; i < m_rank; i++) {
        double center = m_linearSystemWithErrorResult(i,0);
        double width = m_linearSystemWithErrorResult(i,1);
        m_boundingInterval.lowerBounds[i] = std::max(m_boundingInterval.lowerBounds[i], center - width);
        m_boundingInterval.upperBounds[i] = std::min(m_boundingInterval.upperBounds[i], center + width);
        if(m_boundingInterval.lowerBounds[i] >= m_boundingInterval.upperBounds[i]) {
            return 0.0; //Throw out the interval!
        }
    }
    
    //Return the resulting area
    return m_boundingInterval.getArea();
}

template <int Rank>
void IntervalBounder<Rank>::preconditionPolynomials(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations) {
    //Create the polynomials matrix.
    size_t maxPartialSideLength = 0;
    for (size_t i = 0 ; i < _chebyshevApproximations.size(); i++) {
        maxPartialSideLength = std::min(maxPartialSideLength, _chebyshevApproximations[i].getPartialSideLength());
    }
    Eigen::MatrixXd polynomials(m_rank, power(m_rank, maxPartialSideLength));
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
            bool firstPass = true;
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
                
                if(firstPass && spotToInc != 0) {
                    spotToInc = 0;
                }
                firstPass = false;
            }
            inputSpot[spotToInc] = 0;
            spotToInc++;
        }
    }
    
    //Get the preconditioned polynomials
    m_preconditionedErrors = m_linearTermsQR.solve(m_errorTerms); //TOOD: What if this isn't full rank?
    m_preconditionedPolys = m_linearTermsQR.solve(polynomials); //TOOD: What if this isn't full rank?
}

template <int Rank>
double IntervalBounder<Rank>::computeLipshitzConstant(Eigen::VectorXd poly, size_t dim) {
    //poly is of length m_preconditionPolysDegree^m_rank
    //Sum the squares of the spot times the place mod (m_preconditionPolysDegree)^(dim+1) or something like that.
    
}


template <int Rank>
double IntervalBounder<Rank>::updateBoundingIntervalLipshitzSolve(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations) {
    //Get the Lipshitz Constants of the functions.
    Eigen::VectorXd lipshitzConstants;
    for(size_t i = 0; i < m_rank; i++) {
        //lipshitzConstants(i) = computeLipshitzConstant(m_preconditionedPolys.row(i), i); //TODO: Write this
    }
    
    static const double MIN_MOVE = 1e-3;
    bool changed = true;
    while(changed) {
        changed = false;
        for(size_t dim = 0; dim < m_rank; dim++) {
            //Evaluate the dim m_preconditionedPolys at m_boundingInterval.lowerBounds[dim] and m_boundingInterval.upperBounds[dim]
            
            
            //Call get Bound Increase on each function
        }
    }
    
    //TOOD: Finish converting this
    /*#Find M_x and M_y
    LipshitzBounds = []
    for M in Ms:
        LipshitzBounds.append(getLipshitzConstantsND(M))
        
    #Get the mask for the result intervals down each dimension
    intervalMasks = []
    for i in range(len(Ms)):
        mask = [True]*len(Ms)
        mask[i] = False
        intervalMasks.append(mask)
        
    #Iterate until nothing changes. Probably not the fastest way to do it.
    changed = True
    minStep = 1e-3
    while changed:
        changed = False
        #Run it on each poly
        for M, LBounds, error, dimToRun in zip(Ms, LipshitzBounds, errors, np.arange(len(Ms))):
            #Get the reduced dimension Chebyshev Polynomials
            M_hat1 = getLowerDCheb(M, dimToRun, resultIntervals[dimToRun][0])
            M_hat2 = getLowerDCheb(M, dimToRun, resultIntervals[dimToRun][1])

            #Get the other vars down a dimension
            intervalToRun = resultIntervals[intervalMasks[dimToRun]]
            lipshitzBound = LBounds[dimToRun]

            #Optimize over each line
            delta1 = getBoundIncreaseND(M_hat1, error, lipshitzBound, intervalToRun, approxToUse)
            if delta1 > minStep:
                changed = True
                resultIntervals[dimToRun][0] += delta1

            delta2 = getBoundIncreaseND(M_hat2, error, lipshitzBound, intervalToRun, approxToUse)
            if delta2 > minStep:
                changed = True
                resultIntervals[dimToRun][1] -= delta2

            #Check if we've ruled out everything
            if resultIntervals[dimToRun][1] < resultIntervals[dimToRun][0]:
                resultIntervals[:,:] = 0
                return resultIntervals
        
    return resultIntervals*/


    
    return m_boundingInterval.getArea();
}

template <int Rank>
double IntervalBounder<Rank>::computeBoundingInterval(std::vector<ChebyshevApproximation<Rank> >& _chebyshevApproximations) {
    m_boundingInterval.areaFound = false;
    if(updateBoundingIntervalLinearErrorSolve(_chebyshevApproximations) == 0.0) {
        return 0.0;
    }
    return m_boundingInterval.getArea(); //TODO: Remove this once the line check is in.
    m_boundingInterval.areaFound = false;
    
    //TODO: Don't call the constant term check, as the updateBoundingIntervalLinearErrorSolve is better
    //Figure out how to do the absValueSum better.

    
    //Call the Preconditioner.
    preconditionPolynomials(_chebyshevApproximations);
    //Call the Lipshitz Solve
    
    //TODO: Should I just return without the harder checks if it's already be narrowed down to a small enough interval?
    return updateBoundingIntervalLipshitzSolve(_chebyshevApproximations);
}

#endif /* IntervalBounderND_h */
