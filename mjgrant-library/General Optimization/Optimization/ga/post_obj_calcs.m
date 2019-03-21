function [od] = post_fitness_calcs(in,od)
%
% [od] = post_fitness_calcs(in,od)
%
% This function is used to perform any calculations associated with the
% current iteration that require knowledge of the fitness of the individuals.
% This function is executed before the termination criteria is checked.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during 
%                optimization
%
% Output:
%   od - structure containing variables of interest during 
%                optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Elite Individual %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[od.best_fit(1,od.iter),od.best_member(1,od.iter)] = min(od.fit(1,:,od.iter));

return

