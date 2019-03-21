function [in] = set_input
%
% [in] = set_input
%
% This function acts as the interface between the user and the optimization
% algorithm. This function serves as an input deck for all inputs that
% govern the optimization algorithms.
%
% Input:
%   None.
%
% Output:
%   in - optimization algorithm input information [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%
%% Input Parameters %%
%%%%%%%%%%%%%%%%%%%%%%

% Optimization library path
% in.path = 'D:\Documents and Settings\mjgrant.JSC\My Documents\Mike\Scripts\Optimization\Library';
%in.path = 'F:\Scripts\Optimization\Library';
in.path = '.'; % Home

% Determine method
in.method = 'mopso';

% Number of objectives in optimization
in.num_obj = 2;

% Objective information
in.obj_label = {'Obj 1','Obj 2'};

% Determine optimality
in.maximization = [true false];

% Number of individuals in population
in.num_ind = 30;

% Maximum number of iterations to perform
in.max_iter = 100;

% Trade Space Information
in.ts_label = {'Dim 1', 'Dim 2'};

% Upper/lower bounds of design space
in.pos_max = [ 5  5]';
in.pos_min = [-5 -5]';

% Lower bound of velocity
in.vel_min = -1/3*(in.pos_max-in.pos_min);

% Upper bound of velocity
in.vel_max =  1/3*(in.pos_max-in.pos_min);

% Random number generator seed (empty to not specify seed)
in.seed = 0;

% Number of Pareto solutions in external archive
in.num_arch = 50;

% Crowdedness percent - archive portion used for social parameter
in.cp = 0.05;
in.cp = max([in.cp 4.1/in.num_arch]); % use 4 or more points

% Individuality parameter
%in.c1 = 2.0;
in.c1 = 0.0; % GFM: zero'd this, since it doesn't help much for MOPSO

% Sociality parameter
in.c2 = 1.0;

% Constriction factor to limit velocity. Velocity is explicitly limited in code
%  via delta_min and delta_max. Hence, X is not necessary but is retained for
%  future problems. Setting X = 1 does not influence the velocity of the
%  particles.
in.X = 1;

% Inertia Type - Governs type of inertia decay.
% Options: 'cov' = coefficient of variance step decay (need all listed)
%          'rand' = randomly between values (in.ib = [lower_bound upper_bound])
%          'arch' = reduce inertia when archive full
%          'md' = reduce inertia when migration distance less than tolerance
in.w_type = 'md';
in.arch_del_iter = 1;
in.dist = 2;

% Inertia geometric factor (<1)
in.w_factor = 0.975;

% Initial inertia value
in.w_init = 1.0;

% Termination criteria: Percent crowdedness distance -- Determines how much 
% crowdedness distance may vary over the crowdedness distance integral
in.per_cd = 0.02;

% Termination criteria: Determines number of iterations taken in crowdedness
% distance integral
in.del_iter = 20;

% Save frequency (iterations) (set to 0 to not perform intermediate saving)
in.freq_save = in.max_iter*2;

% % Constraint information
% % Name / equality / value
% in.constr = [{'Constraint 1'} {'<='} {0}; ...
%              {'Constraint 2'} {'<='} {0}];

%%%%%%%%%%%%%%%%%%%%
%% Path Execution %%
%%%%%%%%%%%%%%%%%%%%

% Add paths to optimization directory and corresponding algorithm
if ~isdeployed
  addpath(in.path);
  addpath([in.path,'/',in.method]);
end

return

