function [in] = set_input
%
% [in] = set_input
%
% This function acts as the interface between the user and the optimization
% algorithms. This function serves as an input deck for all inputs that
% govern the optimization algorithms.
%
% Input:
%   None.
%
% Output:
%   in - structure containing inputs that are passed to the
%         optimization algorithm
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%
%% Input Parameters %%
%%%%%%%%%%%%%%%%%%%%%%

% Optimization library path
% in.path = 'D:\Documents and Settings\mjgrant.JSC\My Documents\Mike\Scripts\Optimization\Library';
%in.path = 'F:/Work/Scripts/Optimization/Library'; % Flash Drive - Home
in.path = 'C:/Documents and Settings/Mike/Desktop/Work/Scripts/Optimization/Library'; % Home

in.num_obj = 1; % Try to remove later. Need for analyze_results.m

% Determine method
in.method = 'pso';

% Determine optimality
in.maximization = true;

% Number of individuals in population
in.num_ind = 20;

% Maximum number of iterations to perform
in.max_iter = 300;

% Upper bounds of trade-space
in.pos_max = 20*ones(1,2)';

% Lower bounds of trade-space
in.pos_min = -20*ones(1,2)';

% Random number generator seed (empty to not specify seed)
in.seed = 0;

% Termination criteria: Percent fitness -- Determines how much fitness may vary
% over the fitness integral for a given number of iterations.
in.per_fit = 0.03;

% Termination criteria: Determines number of iterations taken in fitness
% integral
in.del_iter = 10;

% Number of neighbors that will communicate with each individual
in.num_nbr = 7;

% Lower bound of velocity
in.vel_min = [-20 -20]';

% Upper bound of velocity
in.vel_max = [20 20]';

% Individuality parameter
in.c1 = 2.0;

% Sociality parameter
in.c2 = 2.0;

% Constriction factor to limit velocity. Velocity is explicitly limited in code
%  via delta_min and delta_max. Hence, X is not necessary but is retained for
%  future problems. Setting X = 1 does not influence the velocity of the
%  particles.
in.X = 1;

% Inertia Type - Governs type of inertia decay.
% Options: 'linear' = linear decay (only need in.w_init)
%          'quadratic' = quadratic decay (only need in.w_init)
%          'cov' = coefficient of variance step decay (need all listed)
in.w_type = 'cov';

% Percentage of top individuals taken for COV calculation
in.cov_per = 0.4;

% Coefficient of variance tolerance to modify inertia
in.cov_tol = 1;

% Inertia geometric factor (<1)
in.w_factor = 0.975;

% Initial inertia value
in.w_init = 1.0;

% Save frequency (iterations) (set to 0 to not perform intermediate saving)
in.freq_save = 5;

% % Constraint information
% % Name / equality / value
% in.constr = [{'Constraint 1'} {'<='} {0}; ...
%              {'Constraint 2'} {'<='} {0}; ...
%              {'Constraint 3'} {'<='} {0}];

% Trade Space Information
in.ts_label = {'Dim 1', 'Dim 2'};

% Objective information
in.obj_label = {'Obj 1'};
           
%%%%%%%%%%%%%%%%%%%%
%% Path Execution %%
%%%%%%%%%%%%%%%%%%%%

% Add paths to optimization directory and corresponding algorithm
if ~isdeployed
  addpath(in.path);
  addpath([in.path,'/',in.method]);
end

return

