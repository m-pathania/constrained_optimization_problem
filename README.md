# constrained_optimization_problem
MATLAB code to find the optimum point (maxima or minima) of a constrained optimization problem

### Functions
**constrv.m :** returns the constraint violation at given point.

**func.m :** function to be optimized. It can return both function value and penalty function value.

**main.m :** main function. Implements the process of constrain based optimization. Performs plotting and also saved the output.

**Marquart.m :** Implementation of Marquardt's method.

**PenatlyFunc.m :** Implementation fo Penalty function method.

**UniD. m :** Performs unidirectional searches using Newton Raphson Method and Bounding Phase Method.

### Files
**input.txt :** first line of the file is a single number representing the question number to be solved.

**OUTPUT.mat :** MALTAB file containing a cell data structure. First column indicates value of R and second column contains a table that stores data for each iteration of marquadt's method for corresponding value of R.

**Report.docx :** Report containing, Problem Definition, Methods used, solutions obtained and observations.