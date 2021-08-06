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

double getLipshitzConstant1D(const double* _array, size_t _arraySize) {
    double value = 0.0;
    for(size_t i = 1; i < _arraySize; i++) {
        value += std::abs(_array[i])*i*i;
    }
    return value;
}

void getLipshitzConstant2D(const double* _array, size_t _arraySize, size_t _arrayFullSize, double& _valueX, double& _valueY) {
    //TODO: Optimize this function
    _valueX = 0.0;
    _valueY = 0.0;
    
    //Iterate through the array
    for(size_t i = 0; i < _arraySize; i++) {
        for(size_t j = 0; j < _arraySize; j++) {
            double arrayVal = std::abs(_array[i*_arrayFullSize + j]);
            _valueX += arrayVal*i*i;
            _valueY += arrayVal*j*j;
        }
    }
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

void evaluateCheb2D(const double* _array, size_t _arraySize, size_t _arrayFullSize, double _value, size_t _dim, std::vector<double>& _result) {
    //Assumes the array is square and _result is initialized to size at least _arraySize
    //Evaluates in dimension _dim
    //This works in place is _array is the array of _result
    if(_arraySize == 0) {
        //Do Nothing
    }
    else if (_arraySize == 1) {
        _result[0] = _array[0];
    }
    else if (_arraySize == 2) {
        if(_dim == 0) {
            _result[0] =  _array[0] + _value * _array[1];
            _result[1] =  _array[_arrayFullSize] + _value * _array[_arrayFullSize + 1];
        }
        else {
            _result[0] =  _array[0] + _value * _array[_arrayFullSize];
            _result[1] =  _array[1] + _value * _array[_arrayFullSize + 1];
        }
    }
    else {
        if(_dim == 0) {
            double doubleValue = 2*_value;
            for(size_t column = 0; column < _arraySize; column++) {
                size_t offset = column * _arrayFullSize;
                double c0 = _array[offset + _arraySize - 2];
                double c1 = _array[offset + _arraySize - 1];
                //Run the loop
                double temp;
                for(size_t i = _arraySize-3; i < _arraySize; i--) {
                    temp = c0;
                    c0 = _array[offset + i] - c1;
                    c1 = temp + c1 * doubleValue;
                }
                _result[column] = c0 + c1*_value;
            }
        }
        else {
            double doubleValue = 2*_value;
            for(size_t column = 0; column < _arraySize; column++) {
                size_t offset = _arrayFullSize * (_arraySize-1);
                double c1 = _array[column + offset];
                offset -= _arrayFullSize;
                double c0 = _array[column + offset];
                //Run the loop
                double temp;
                while(offset > _arrayFullSize) {
                    offset -= _arrayFullSize;
                    temp = c0;
                    c0 = _array[column + offset] - c1;
                    c1 = temp + c1 * doubleValue;
                }
                _result[column] = c0 + c1*_value;
            }
        }
    }

}










#endif /* BoundingIntervalUtilities_h */
