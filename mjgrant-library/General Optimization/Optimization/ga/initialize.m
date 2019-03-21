function [od] = initialize(in,od)
%
% [od] = initialize(in,od)
%
% This function is used to allocate memory locations for matrices and
% initialize scalar values.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during 
%         optimization
%
% Output:
%   od - structure containing variables of interest during 
%         optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trade Space Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine number of trade space dimensions from bounds
dim1 = length(in.pos_max); % Length of trade space upper bound vector
dim2 = length(in.pos_min); % Length of trade space lower bound vector

% Assign trade space dimension and error check for inconsistent boundary vector
% lengths
if dim1 ~= dim2
  % Upper and lower bounds do not represent the same number of trade space
  % dimensions
  error('Upper and lower bound matrices must be same size.');
else
  % Upper and lower bound vectors same size, assign trade space dimension value
  od.dim = dim1;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Random Number Seed %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if seed assigned to random number generator
if ~isempty(in.seed)
  % Reset random number seed if want to have same initial conditions
  rand('state',in.seed);
end

%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Values %%
%%%%%%%%%%%%%%%%%%%%%%%

% Matrix used to store the fitness history of each individual
od.fit = NaN*ones(1,in.num_ind,in.max_iter);

% Matrix used to store the position history of all individuals
od.pop = NaN*ones(od.dim,in.num_ind,in.max_iter);

% Optimization execution time
od.time = NaN;

% Trade space dimenstions
od.dim = od.dim;

% Vector used to store the best fitness of each iteration
od.best_fit = NaN*ones(1,in.max_iter);

% Matrix of parent selection history
od.parents = NaN*ones(od.dim,in.num_ind,in.max_iter);

% Vector of best member in population history
od.best_member = NaN*ones(1,in.max_iter);

% Constraint matrix
if isfield(in,'constr')
  od.num_constr = length(in.constr(:,1));
else
  od.num_constr = 0;
end
od.constr = NaN*ones(od.num_constr,in.num_ind,in.max_iter);

return

