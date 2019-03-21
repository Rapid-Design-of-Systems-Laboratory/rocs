function [fitness,constr] = fitness_func(mat_input)
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
%   mat_input - Matrix that assumes that each row specifies a dimension
%               that is varied to find the optimal solution. Thus, each 
%               column represents the location of the individual.
%
% Output:
%   fitness - row vector containing all of fitness values for each individual
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%
%% Assign Input %%
%%%%%%%%%%%%%%%%%%

X = mat_input;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Unconstrained Problems %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Bumpy surface
fitness = -X(1,:).^2 - X(2,:).^2 + 2*X(1,:) + X(2,:) + ...
  100*sin(X(1,:)) + 200*sin(X(2,:));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Constrained Optimization Test Cases %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Test Case 1 - Bounds Unknown
% fitness = (X(1,:) - 2).^2 + (X(2,:) - 1).^2;
% constr(1,:) = 2*X(2,:) - 1 - 0.5 - X(1,:);
% constr(2,:) = X(1,:) - (2*X(2,:) - 1 + 0.5);
% constr(3,:) = X(1,:).^2/4 + X(2,:).^2 - 1;

return

