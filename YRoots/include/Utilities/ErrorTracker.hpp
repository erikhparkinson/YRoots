//
//  ErrorTracker.hpp
//  YRoots
//
//  Created by Erik Hales Parkinson on 8/05/21.
//  Copyright © 2021 Erik Hales Parkinson. All rights reserved.
//

#ifndef ErrorTracker_h
#define ErrorTracker_h

constexpr double machineEpsilon = std::numeric_limits<double>::epsilon();

struct ErrorTracker {
    ErrorTracker() : value(0), error(0) {}
    ErrorTracker(double v, double e) : value(v), error(e) {}
    //Assume doubles have error in the last digit, unless it is an exact integer.
    ErrorTracker(double v) : value(v), error(std::abs(v)*machineEpsilon) {
        if(value == static_cast<int64_t>(value)) {
            error = 0;
        }
    }
    //Cast integers with 0 error, as they are exact.
    ErrorTracker(short int v) : value(v), error(0) {}
    ErrorTracker(unsigned short int v) : value(v), error(0) {}
    ErrorTracker(int v) : value(v), error(0) {}
    ErrorTracker(unsigned int v) : value(v), error(0) {}
    ErrorTracker(long int v) : value(v), error(0) {}
    ErrorTracker(unsigned long int v) : value(v), error(0) {}
    ErrorTracker(long long int v) : value(v), error(0) {}
    ErrorTracker(unsigned long long int v) : value(v), error(0) {}

    //For printing
    friend std::ostream& operator<<(std::ostream& strm, const ErrorTracker& obj) {
        strm << "Value: " << obj.value << ", Error: " << obj.error;;
        return strm;
    }

    //Override addition
    ErrorTracker& operator+=(const ErrorTracker& rhs)
    {
        error = std::max(std::max(std::abs(value), std::abs(rhs.value))* machineEpsilon, std::max(error, rhs.error));
        value += rhs.value;
        return *this; // return the result by reference
    }
    friend ErrorTracker operator+(ErrorTracker lhs, const ErrorTracker& rhs)
    {
        lhs += rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }

    //Override subtraction
    ErrorTracker& operator-=(const ErrorTracker& rhs)
    {
        error = std::max(std::max(std::abs(value), std::abs(rhs.value))* machineEpsilon, std::max(error, rhs.error));
        value -= rhs.value;
        return *this; // return the result by reference
    }
    friend ErrorTracker operator-(ErrorTracker lhs, const ErrorTracker& rhs)
    {
        lhs -= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }

    //Override multiplication
    ErrorTracker& operator*=(const ErrorTracker& rhs)
    {
        error = std::abs(value) * rhs.error + error * std::abs(rhs.value) + error * rhs.error;
        value *= rhs.value;
        return *this; // return the result by reference
    }
    friend ErrorTracker operator*(ErrorTracker lhs, const ErrorTracker& rhs)
    {
        lhs *= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }

    //Override Division
    ErrorTracker& operator/=(const ErrorTracker& rhs)
    {
        if (rhs.error >= std::abs(rhs.value)) {
            throw std::runtime_error("Divided by interval 0! Function can't be evaulated stabily.");
        }
        const double val1 = (value + error) / (rhs.value + rhs.error); //TODO: Improve this!
        const double val2 = (value + error) / (rhs.value - rhs.error);
        const double val3 = (value - error) / (rhs.value + rhs.error);
        const double val4 = (value - error) / (rhs.value - rhs.error);
        value /= rhs.value;
        error = std::max(std::max(std::abs(val1-value), std::abs(val2 - value)), std::max(std::abs(val3 - value), std::abs(val4 - value)));
        return *this; // return the result by reference
    }
    friend ErrorTracker operator/(ErrorTracker lhs, const ErrorTracker& rhs)
    {
        lhs /= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }

    double value;
    double error;
};

//Evaulations of functions with error tracking. For some input x with error e, the error of f(x) will be bounded by abs(f(x))*machineEpsilon + lipshitzConstant * e. The lipshitzConstant only needs to apply on the domain [x-e.x+e], so for small e, f'(x) will give a good approximation of this.
ErrorTracker sin(ErrorTracker x) {
    const double lipshitzConstant = std::abs(cos(x.value)); //TODO: Make a more rigorous Lipshitz Constant
    const double eval = sin(x.value);
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker cos(ErrorTracker x) {
    const double lipshitzConstant = std::abs(sin(x.value)); //TODO: Make a more rigorous Lipshitz Constant
    const double eval = cos(x.value);
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker tan(ErrorTracker x) {
    return sin(x) / cos(x);
}

ErrorTracker sinh(ErrorTracker x) {
    const double lipshitzConstant = std::abs(cosh(x.value)); //TODO: Make a more rigorous Lipshitz Constant
    const double eval = sinh(x.value);
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker cosh(ErrorTracker x) {
    const double lipshitzConstant = std::abs(sinh(x.value)); //TODO: Make a more rigorous Lipshitz Constant
    const double eval = cosh(x.value);
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker tanh(ErrorTracker x) {
    return sinh(x) / cosh(x);
}

ErrorTracker exp(ErrorTracker x) {
    const double eval = exp(x.value);
    const double lipshitzConstant = eval; //TODO: Make a more rigorous Lipshitz Constant
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker sqrt(ErrorTracker x) {
    const double eval = sqrt(x.value);
    const double lipshitzConstant = 0.5/eval; //TODO: Make a more rigorous Lipshitz Constant
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker log(ErrorTracker x) {
    const double eval = log(x.value);
    const double lipshitzConstant = 1 / x.value; //TODO: Make a more rigorous Lipshitz Constant
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker log2(ErrorTracker x) {
    const double eval = log2(x.value);
    const double lipshitzConstant = log(2) / x.value; //TODO: Make a more rigorous Lipshitz Constant
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker log10(ErrorTracker x) {
    const double eval = log10(x.value);
    const double lipshitzConstant = log(10) / x.value; //TODO: Make a more rigorous Lipshitz Constant
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstant);
}

ErrorTracker pow(ErrorTracker x, size_t y) {
    const double eval = pow(x.value, y);
    const double lipshitzConstantX = std::abs(y * eval / x.value); //TODO: Make a more rigorous Lipshitz Constant
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstantX);
}

ErrorTracker pow(ErrorTracker x, ErrorTracker y) {
    if (x.value < 0 && y.error != 0) {
        throw std::runtime_error("Taking the power of a negative can't be evaluated stabily!");
    }

    const double eval = pow(x.value, y.value);
    const double lipshitzConstantX = std::abs(y.value * eval / x.value); //TODO: Make a more rigorous Lipshitz Constant
    const double lipshitzConstantY = std::abs(log(x.value) * eval); //TODO: Make a more rigorous Lipshitz Constant
    return ErrorTracker(eval, std::abs(eval)*machineEpsilon + x.error * lipshitzConstantX + (y.error == 0 ? 0 : y.error * lipshitzConstantY));
}

#endif //ErrorTracker_h
