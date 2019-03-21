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
%            optimization algorithm
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%% WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%% GA CURRENTLY NOT SET UP TO PERFORMED CONSTRAINED OPTIMIZATION. ISSUE OCCURS
%%% IN ADVANCE_POP.M DURING PARENT SELECTION!!!! CAN'T SELECT NAN PARENTS!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%
%% Input Parameters %%
%%%%%%%%%%%%%%%%%%%%%%

% Optimization library path
% in.path = 'D:\Documents and Settings\mjgrant.JSC\My Documents\Mike\Scripts\Optimization\Library';
% in.path = 'G:\Work\Scripts\Optimization\Library';
in.path = 'C:/Documents and Settings/Mike/Desktop/Work/Scripts/Optimization/Library'; % Home

in.num_obj = 1; % Try to remove later. Need for analyze_results.m

% Determine method
in.method = 'ga';

% Determine optimality
in.maximization = true;

% Number of individuals in population
in.num_ind = 20;

% Maximum number of iterations to perform
in.max_iter = 300;

% Upper bounds of trade-space (inf valid)
in.pos_max = [20 20]';

% Lower bounds of trade-space (-inf valid)
in.pos_min = [-20 -20]';

% Random number generator seed (empty to not specify seed)
in.seed = 0;

% Termination criteria: Percent fitness -- Determines how much fitness may vary
% over the fitness integral for a given number of iterations.
in.per_fit = 0.02;

% Termination criteria: Determines number of iterations taken in fitness
% integral
in.del_iter = 30;

% Elitism on/off
in.elitism = true;

% Probability of crossover for each mating pairs of individuals
in.cross_prob = 0.80;

% Probability of mutation for each individual
in.mutate_prob_ind = 1.0;

% Probability of mutation for each trait
in.mutate_prob_trait = 0.05;

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

