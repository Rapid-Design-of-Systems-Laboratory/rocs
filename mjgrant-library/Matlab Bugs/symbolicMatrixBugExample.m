function symbolicMatrixBugExample

% Both of the following blocks work in R2008b, but the first no longer
% works with R2010a.

clc

syms a

% Error with symbolic nan
val = sym(nan(1,1))
val(1,1) = a % This line doesn't change val in R2010a but does in R2008b
val = a % This line does change val in both R2010a and R2008b

% Works with symbolic zeros
val = sym(zeros(1,1))
val(1,1) = a % This line does change val in both R2010a and R2008b

return

