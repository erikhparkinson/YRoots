# YRoots
Implementing BYU's YRoots Rootfinding code in C++

# Running
Download the codebase as a zip file, then open a command line and in the main folder run `make`.
This will create the executable `yroots-solver`.
Yroots can be run from the command line as `./yroots <input file>`. If no input file is given `input.txt` is assumed.

# File Input

The input file is of the following format. All whitespace is ignored by the parser.
Any line that starts with a `#` is ignored by the parser.

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
* targetTol : Defaults to 1e-15
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
[-pi, cos(sin(ln(7)))];
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

This must exists for every function declared in the first line, but can also exist for other subfunctions that are then used in the main funcitons. 

The function parsing was set up to hopefully be fairly intuitive, and you should be able to just write function as you'd think it should be. The typical order of operations if followed.
Specifically, the function parser works recursively with the following supported function formats.
* `<subfunction1> + <subfunction2>` - Addition. Example: `x+y`.
* `<subfunction1> - <subfunction2>` - Subtraction. Example: `x-y`.
* `<subfunction1> * <subfunction2>` - Multiplication. Example: `x*y`.
* `<subfunction1> / <subfunction2>` - Division. Example: `x/y`.
* `<subfunction1> ^ <subfunction2>` - Power. The Number is optional. Example: `3x^2`.
* `<subfunction1> ** <subfunction2>` - Power. Same as above but makes copying functions from Python easier.  Example: `3x**2`.
* `<Expression>(<subfunction>)` - More complicated expressions as defined below.
* `<Variable Name>` - The name of one of the variable inputs.
* `<Function Name>` - The name of one of the funtion inputs. This function must have been previously defined.
* `<Number>` - Numbers.
* `<Constants>` - Numeric constants. Currently `e` and `pi` are supported, not case sensitive.
* `(<subfunction>)` - Using parenthesis.

The supported `Expression` for complex functions with one input are
* `sin`, `cos`, `tan`,`sinh`, `cosh`, `tanh` - Trigonimetric Functions. Example: `3sin(x)`.
* `sqrt`,`exp` - Square root, Exponential Function, and Natural Logarithm.  Example: `-2.5exp(x)`.
* `log`,`log2`,`log10` - Logarithm base e, 2, 10. Example: `3log(x)`.
* `prod`,`sum` - Product and sum notation. `sum(<func>,<var>,start,end)` will sum `func` for the values of `[start,end]` for `var`. The syntax for `prod` is the same, but it will take the product of the functions. `start` and `end` must be integerts, it will be included in the sum/product. Examples: `prod(sin(x+i),i,0,10)`, `sum(x^i,i,-50,-1).`
* `T<n>` - Chebyshev Polynomials. gives the nth Chebyshev polynomials of the first kind. Example: `T2(x)`.
 
# File Output

## Roots
The roots will be outputted into a file called roots.csv. Each line will contain one root, with the coordinate in each dimension seperated by a comma.

## Intervals
If the option trackIntervals is set true, the intervals that were solved will be printed to the file intervals.txt.

## Timing
If the option useTimer is set true, timing information about the solve will be printed to the file timing.txt.
