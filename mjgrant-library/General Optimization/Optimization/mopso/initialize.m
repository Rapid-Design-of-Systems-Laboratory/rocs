function [od] = initialize(in,od)
%
% [od] = initialize(inputs,od)
%
% This function is used to allocate memory locations for matrices and
% initialize scalar values.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - structure containing variables of interest during 
%                optimization
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

% Matrix used to store the objective history of each individual
od.fit = NaN*ones(in.num_obj,in.num_ind,in.max_iter);

% Matrix used to store objective history of Pareto-optimal solutions
od.fit_arch = NaN*ones(in.num_obj,in.num_arch,in.max_iter);

% Matrix used to store position history of Pareto-optimal solutions
od.pos_arch = NaN*ones(od.dim,in.num_arch,in.max_iter);

% Matrix used to store history of position of all individuals
od.pop = NaN*ones(od.dim,in.num_ind,in.max_iter);

% Matrix used to store history of velocity of all individuals
od.vel = NaN*ones(od.dim,in.num_ind,in.max_iter);

% Vector used to store history of cognitive parameter
od.r1 = NaN*ones(1,in.max_iter);

% Vector used to store history of social parameter
od.r2 = NaN*ones(1,in.max_iter);

% Vector used to store history of inertia
od.w = NaN*ones(1,in.max_iter);

% Matrix used to store the location used in the cognitive parameter of the
% velocity equation. Note that "best" does not necessarily exist in a
% multiobjective problem.
od.best = NaN*ones(od.dim,in.num_ind,in.max_iter);

% Matrix used to store the corresponding objective values. Note that "best" 
% does not necessarily exist in a multiobjective problem.
od.best_fit = NaN*ones(in.num_obj,in.num_ind,in.max_iter);

% Optimization execution time
od.time = NaN;

% Trade space dimenstions
od.dim = od.dim;

% Vector used to store the history of crowdedness distance
od.cd = NaN*ones(1,in.num_arch,in.max_iter);

% Vector used to store the history of average crowdedness distance
od.cd_mean = NaN*ones(1,in.max_iter);

% Vector used to store history of coefficient of variance
od.cov = NaN*ones(1,in.max_iter);

% Termination flag
od.flag = NaN*ones(1,in.max_iter); % NaN = everything normal

% Matrix used to save all objective information, regardless if constraints were
% violated.
od.fit_save = NaN*ones(in.num_obj,in.num_ind,in.max_iter);

% Migration distance of Pareto front
if strcmp(in.w_type,'md')
  od.md = NaN*ones(1,in.max_iter);
end

% Determine number of constraints
if isfield(in,'constr')
  od.num_constr = length(in.constr(:,1));
else
  od.num_constr = 0;
end

% Matrix used to store history of constraint information
od.constr = NaN*ones(od.num_constr,in.num_ind,in.max_iter);

return

