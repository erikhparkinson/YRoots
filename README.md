# YRoots
Implementing BYU's YRoots Rootfinding code in C++

# Running
Yroots can be run from the command line with the following syntax.

./yroots \<input file\> \<output file\>

# File Input

The input file is of the following format. All whitespace is ignored by the parser.

```
PARAMETERS;
<Parameter Information>
PARAMETERS END;

INTERVAL;
<Parameter Information>
INTERVAL END;

FUNCTIONS;
<Function Information>
FUNCTIONS END;

END;
```

An example input file is
```
PARAMETERS;
numThreads = 3;
PARAMETERS END;

INTERVAL;
[-1, 1/pi];
[ln(2), e];
INTERVAL END;

FUNCTIONS;
function f1, f2;
variable_group x1, x2;
a = sqrt(x1^2 + x2^2)
f1 = sin(x1) + 2 + ln(a);
f2 = cos(a)/(log(x2, 10)) + 2 + ln(x1 + 3);
FUNCTIONS END;

END;
```

## Parameters
Each Parameter Line should be of the form
```<name> = <value>;```

Currently supported parameters are
* numThreads
  * Must be an integer. -1 will use the maximum number of possible threads determined by `std::thread::hardware_concurrency()`. Otherwise this value must be positive. Will be automatically be capped at `std::thread::hardware_concurrency()`.
  
 
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
  * `+`, `-`, `*`, `/`, and `^`.
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


