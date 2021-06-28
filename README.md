# YRoots
Implementing BYU's YRoots Rootfinding code in C++

# Running
Yroots can be run from the command line with the following syntax.

./yroots \<input file\>

# File Input

The input file is of the following format. All whitespace is ignored by the parser.

```
PARAMETERS;
<Parameter Information>
PARAMETERS_END;

INTERVAL;
<Parameter Information>
INTERVAL_END;

FUNCTIONS;
<Function Information>
FUNCTIONS_END;

END;
```

An example input file is
```
PARAMETERS;
numThreads = 3;
PARAMETERS_END;

INTERVAL;
[-1, 1/pi];
[ln(2), e];
INTERVAL_END;

FUNCTIONS;
function f1, f2;
variable_group x1, x2;
a = sqrt(x1^2 + x2^2)
f1 = sin(x1) + 2 + ln(a);
f2 = cos(a)/(log(x2, 10)) + 2 + ln(x1 + 3);
FUNCTIONS_END;

END;
```

## Parameters
This section is optional. If not used everything will use the default values. Every field is also optional. 
Each Parameter Line should be of the form
```<name> = <value>;```

Currently supported parameters are
* numThreads : Defaults to 1
  * Must be an integer. Defaults to 1. -1 will use the maximum number of possible threads determined by `std::thread::hardware_concurrency()`. Otherwise this value must be positive.
* relApproxTol : Defaults to 1e-10
  * Must be a positive number. Used with absApproxTol to determine if an approximation is good.
* absApproxTol : Defaults to 1e-10
  * Must be a positive number. An approximation if good if `error < absApproxTol + inf_norm * relApproxTol`.
* goodZerosFactor : Defaults to 100
  * Must be a non-negative number. A zero is declared to be on the real line [-1,1] if it's real and imaginary parts are no more than tol * approx_error away from it.
* minGoodZerosTol: Defaults to 1e-5
  * The min value of tol * approx_error in the goodZerosFactor calculation.
* approximationDegree : Defaults to 20
  * The initial approximation degree used.
* maxLevel : Defaults to 50
  * How many levls deep the subdivision will go before giving up.
* trackIntervals : Defaults to true
  * If true the intervals are tracked and the results of how each one is solved is saved in `intervals.txt`.
* useTimer : Defaults to false
  * If true timing details are recorded and saved in `timing.txt`.
 
## Intervals
This defines the intervals on which the roots are found. Currently this is only defined for an n-dimensional box. This section should contain a line for each of the _n_ variables defined in the _Functions_ section. Each line is of the form
```[<lowerBound>, <upperBound>];```

The bounds can be any type of constant numberical expression that can be interpreted by _Functions_ as defined in the next section.
For example, all of the following are valid intervals.
```
[-11, 1.5];
[2^(-3), e];
[-pi, cos(sin(log(7)))];
```

## Functions
The first line is a declaration of the function names. It is of the form
```
function <name_1>, <name_2>, ... <name_n>;
```

The second line is a declaration of the variable names. It is of the form
```
variable_group <var_1>, <var_2>, ... <var_n>;
```

The reset of the lines all declare the function. They are of the form
```
<function_name> = <expression>;
```

This must exists for every function declared in the first line, but can also exist for other subfunctions that are then used in the main funcitons. Functions can contain the following types. Letters are not cases sensitive. 
* Constants
  * These include numbers, `e`, and `pi`.
* Addition, Subtraction, Multiplication, Division and Powers
  * `+`, `-`, `*`, `/`, and `^`. To make Python users happy, `**` is interpreted as `^`.
* Trig functions
  * `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`. Used as `sin(<expression>)`.
* Square roots
  * `sqrt(<expression>)`.
* Exponentials and Logirithms
  * `exp(<expression>)`, `log(<expression>, <b>)`, and `ln(<expression>)`. ln is the base e logarithm. log is the base b logarithm, where b is not neccesarily constant.
* Chebyshev Polynomials
  * `T<n>(<expression>)` gives the nth Chebyshev polynomials of the first kind of \<expression\>, 
* Variables and Other Functions
  * The name of a defined variable or function. Any function is assumed to be of the same dimension and take in the same variables.


