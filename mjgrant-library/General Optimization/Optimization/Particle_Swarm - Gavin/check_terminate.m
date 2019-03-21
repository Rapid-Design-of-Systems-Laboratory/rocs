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
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%   terminate - boolean signifying if the optimization algorithm should 
%               terminate
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Terminate boolean
terminate = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine if Pareto Front Converging %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if Pareto front converging. The average crowdedness distance is 
% integrated over the number of iterations specified in the input.  If the 
% integrated average crowdedness distance is less than the percentage 
% crowdedness distance specified in the input, the optimization process is 
% terminated.
if (od.iter > in.del_iter)
  
  % Check integral of average crowdedness distance
  integral = sum(abs(diff((od.cd_mean(od.iter-in.del_iter:od.iter)))));
  
  % If integral less than tolerance specified by user, terminate optimization
  if integral <= in.per_cd*od.cd_mean(od.iter)
    terminate = true;
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Protect Premature Termination %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If termination is commanded (but archive is not full) then the archive has
% probably not grown in the number of iterations used for the convergence
% criteria. Force inertia reduction.

% Determine size of archive
[R,C] = find(isnan(od.fit_arch(:,:,od.iter)) == 0);
C = unique(C);

if (terminate) && (length(C) ~= in.num_arch)
  % Premature convergence, archive not full
  od.flag(od.iter) = 1;
  terminate = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Maximum Number of Iterations Criteria %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if optimization algorithm reached maximum number of iterations
if (od.iter == in.max_iter)
  terminate = true;
end

return
