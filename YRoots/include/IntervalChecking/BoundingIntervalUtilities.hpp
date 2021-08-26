//
//  BoundingIntervalUtilities.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 12/23/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef BoundingIntervalUtilities_h
#define BoundingIntervalUtilities_h

#include <stddef.h>
#include <cmath>
#include <vector>
#include <math.h>
#include "Utilities/utilities.hpp"

//Helpers for the 1-D check

double getLipshitzConstant1D(const double* _array, size_t _arraySize) {
    double value = 0.0;
    for(size_t i = 1; i < _arraySize; i++) {
        value += std::abs(_array[i])*i*i;
    }
    return value;
}

double evaluateCheb1D(const double* _array, size_t _arraySize, double _value) {
    if(_arraySize == 0) {
        return 0;
    }
    else if (_arraySize == 1) {
        return _array[0];
    }
    else if (_arraySize == 2) {
        return _array[0] + _value * _array[1];
    }
    else {
        double doubleValue = 2*_value;
        double c0 = _array[_arraySize - 2];
        double c1 = _array[_arraySize - 1];
        double temp;
        for(size_t i = _arraySize-3; i < _arraySize; i--) {
            temp = c0;
            c0 = _array[i] - c1;
            c1 = temp + c1 * doubleValue;
        }
        return c0 + c1*_value;
    }
}

//Helpers for the 2-D check

double getOptimizeLineReturnValue(double minVal, double maxVal, double error) {
    if (minVal > error) {
        return minVal - error;
    }
    else if (maxVal < -error) {
        return -maxVal - error;
    }
    else {
        return 0;
    }
}

double optimizeLine0(const Eigen::VectorXd& poly, double low, double high) {
    double value = poly(0);
    double error = 0.0;
    for(int i = 1; i < poly.size(); i++) {
        error += std::abs(poly(i));
    }
    return getOptimizeLineReturnValue(value, value, error);
}

double optimizeLine1(const Eigen::VectorXd& poly, double low, double high) {
    if(unlikely(poly(1) == 0.0)) {
        return optimizeLine0(poly, low, high);
    }
    
    double b = poly(0);
    double a = poly(1);
    double error = 0.0;
    for(int i = 2; i < poly.size(); i++) {
        error += std::abs(poly(i));
    }
    double v1 = b + low * a;
    double v2 = b + high * a;
    double minValue = std::min(v1, v2);
    double maxValue = std::max(v1, v2);
    return getOptimizeLineReturnValue(minValue, maxValue, error);
}

double inline chebValQuad(double c, double b, double a, double value) {
    return c + (b * value) + (a * (2 * value * value - 1));
}

double optimizeLine2(const Eigen::VectorXd& poly, double low, double high) {
    if(unlikely(poly(2) == 0.0)) {
        return optimizeLine1(poly, low, high);
    }

    double c = poly(0);
    double b = poly(1);
    double a = poly(2);
    
    double error = 0.0;
    for(int i = 3; i < poly.size(); i++) {
        error += std::abs(poly(i));
    }
    double v1 = chebValQuad(c, b, a, low);
    double v2 = chebValQuad(c, b, a, high);
    double minValue = std::min(v1, v2);
    double maxValue = std::max(v1, v2);
    if(a != 0) {
        double der0point = -b/(4*a);
        if (der0point > low && der0point < high) {
            double v = chebValQuad(c, b, a, der0point);
            minValue = std::min(minValue, v);
            maxValue = std::max(maxValue, v);
        }
    }
    return getOptimizeLineReturnValue(minValue, maxValue, error);
}

double inline chebValCubic(double d, double c, double b, double a, double value) {
    const double temp = 2 * value * value;
    return d + (c * value) + (b * (temp - 1)) + (a * (2 * temp * value - 3 * value));
}

double optimizeLine3(const Eigen::VectorXd& poly, double low, double high) {
    if(unlikely(poly(3) == 0.0)) {
        return optimizeLine2(poly, low, high);
    }

    double d = poly(0);
    double c = poly(1);
    double b = poly(2);
    double a = poly(3);
    double error = 0.0;
    for(int i = 4; i < poly.size(); i++) {
        error += std::abs(poly(i));
    }
    double v1 = chebValCubic(d, c, b, a, low);
    double v2 = chebValCubic(d, c, b, a, high);
    double minValue = std::min(v1, v2);
    double maxValue = std::max(v1, v2);
    
    if(unlikely(b == 0.0)) { //Solving 12ax^2-3a+c = 0
        double value = 0.25 - c/(12*a);
        if(value > 0) {
            double der0Point = sqrt(value);
            double v = chebValCubic(d, c, b, a, der0Point);
            minValue = std::min(minValue, v);
            maxValue = std::max(maxValue, v);
            v = chebValCubic(d, c, b, a, -der0Point);
            minValue = std::min(minValue, v);
            maxValue = std::max(maxValue, v);
        }
    }
    else { //Solving 12ax^2 + 4bx -3a+c = 0
        double discr = b*b - 3*a*(c-3*a);
        if (discr >= 0) {
            //When solving ax^2 + bx + c, discr = b^2 - 4ac
            //Solution 1 at (-b - sign(b) * sqrt(discr)) / (2a)
            //Solution 2 at (2c) / (-b - sign(b) * sqrt(discr))
            double temp = -b + (b > 0 ? -sqrt(discr) : sqrt(discr));
            double der0point1 = temp / (6*a);
            double der0point2 = (c - 3*a) / (2*temp);

            if(der0point1 > low && der0point1 < high) {
                double v = chebValCubic(d, c, b, a, der0point1);
                minValue = std::min(minValue, v);
                maxValue = std::max(maxValue, v);
            }
            if(der0point2 > low && der0point2 < high) {
                double v = chebValCubic(d, c, b, a, der0point2);
                minValue = std::min(minValue, v);
                maxValue = std::max(maxValue, v);
            }
        }
    }
    
    return getOptimizeLineReturnValue(minValue, maxValue, error);
}




#endif /* BoundingIntervalUtilities_h */
