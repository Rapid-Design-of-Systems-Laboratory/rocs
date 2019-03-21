function [fitness,constr] = eval_pop(mat_input)
%
% [fitness] = fitness_func(mat_input)
%
% This function is used to return the fitness of all of the test solutions.
% A simulation could be performed in this function or a continuous function
% can be tested. The output returned should represent a fitness parameter
% that will be maximized. In order to perform a minimization process,
% negate the fitness values.
%
% Input:
%   mat_input - 2-D matrix that assumes that each row specifies a dimension
%               that is varied to find the optimal solution. Thus, each 
%               column represents the location of the individual.
%
% Output:
%   fitness - row vector containing all of fitness values for each individual
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%
%% Assign input %%
%%%%%%%%%%%%%%%%%%

X = mat_input;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unconstrained Problems %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bumpy surface
fitness = -X(1,:).^2 - X(2,:).^2 + 2*X(1,:) + X(2,:) + ...
  100*sin(X(1,:)) + 200*sin(X(2,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained Optimization Test Cases %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Test Case 1 - Bounds Unknown
% fitness = (X(1,:) - 2).^2 + (X(2,:) - 1).^2;
% constr(1,:) = 2*X(2,:) - 1 - 0.5 - X(1,:);
% constr(2,:) = X(1,:) - (2*X(2,:) - 1 + 0.5);
% constr(3,:) = X(1,:).^2/4 + X(2,:).^2 - 1;

% I1 = find(X < 2*Y - 1 - 0.5);
% I2 = find(X > 2*Y - 1 + 0.5);
% I3 = find(X.^2/4 + Y.^2 - 1 > 0);
% I_sub = union(I1,I2);
% I = union(I_sub,I3);
% fitness(I) = NaN;

% %% Test Case 2 - 
% fitness = (X - 10).^3 + (Y - 20).^3;
% I1 = find(100 - (X - 5).^2 - (Y - 5).^2 > 0);
% I2 = find((X - 6).^2 + (Y - 5).^2 - 82.81 > 0);
% I = union(I1,I2);
% fitness(I) = NaN;

% %% Test Case 3 - -10 <= xi <= 10
% fitness = 5*(X1 + X2 + X3 + X4) - 5*(X1.^2 + X2.^2 + X3.^2 + X4.^2) - ...
%   (X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13);
% I1 = find(2*X1 + 2*X2 + X10 + X11 > 10);
% I2 = find(2*X1 + 2*X3 + X10 + X12 > 10);
% I3 = find(2*X2 + 2*X3 + X11 + X12 > 10);
% I4 = find(-8*X1 + X10 > 0);
% I5 = find(-8*X2 + X11 > 0);
% I6 = find(-8*X3 + X12 > 0);
% I7 = find(-2*X4 - X5 + X10 > 0);
% I8 = find(-2*X6 - X7 + X11 > 0);
% I9 = find(-2*X8 - X9 + X12 > 0);
% 
% I = union(I1,I2);
% I = union(I,I3);
% I = union(I,I4);
% I = union(I,I5);
% I = union(I,I6);
% I = union(I,I7);
% I = union(I,I8);
% I = union(I,I9);
% 
% fitness(I) = NaN;
% 
% %% Test Case 4
% fitness = (X1 - 10).^2 + 5*(X2 - 12).^2 + X3.^4 + 3*(X4 - 11).^2 + ...
%            10*X5.^6 + 7*X6.^2 + X7.^4 - 4*X6.*X7 - 10*X6 - 8*X7;
% I1 = find(127 - 2*X1.^2 - 3*X2.^4 - X3 - 4*X4.^2 - 5*X5 < 0);
% I2 = find(282 - 7*X1 - 3*X2 - 10*X3.^2 - X4 + X5 < 0);
% I3 = find(196 - 23*X1 - X2.^2 - 6*X6.^2 + 8*X7 < 0);
% I4 = find(-4*X1.^2 - X2.^2 + 3*X1.*X2 - 2*X3.^2 - 5*X6 + 11*X7 < 0);
% I = union(I1,I2);
% I = union(I,I3);
% I = union(I,I4);
% fitness(I) = NaN;


return

