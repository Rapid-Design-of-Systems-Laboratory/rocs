function [od,terminate] = check_terminate(in,od)
%
% [terminate] = check_terminate(in,od)
%
% This function determines if the optimization process should be
% terminated. Both the user inputs and optimization data are passed to this
% function in order to provide various means by which to terminate the
% optimization algorithm.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during 
%         optimization
%
% Output:
%   terminate - boolean signifying if the optimization algorithm should
%                terminate
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Terminate boolean
terminate = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitness Convergence Criteria %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if fitness converging. The fitness is integrated over the number of
% iterations specified in the input.  If the integrated fitness is less than the
% percentage fitness specified in the input, the optimization process is 
% terminated.  Use <= so if optimal fitness value is zero, might converge.
if (od.iter > in.del_iter) && (sum(abs(diff((od.best_fit(od.iter-...
    in.del_iter:od.iter))))) <= in.per_fit*abs(od.best_fit(od.iter)))
  terminate = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Maximum Number of Iterations Criteria %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if optimization algorithm reached maximum number of iterations
if (od.iter == in.max_iter)
  terminate = true;
end

return