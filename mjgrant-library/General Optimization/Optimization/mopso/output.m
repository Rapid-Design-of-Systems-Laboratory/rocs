function [] = output(in,od,ui)
%
% [] = output(in,od)
%
% This function is used to provide output immediately after the multi-objective
% optimization process has completed.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   None. Although no values are returned to a calling function, output is
%   printed to the screen and the results are stored in results.mat.
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Feb. 2006


%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitness Conversion %%
%%%%%%%%%%%%%%%%%%%%%%%%

% If performing maximization, objective values were negated. Need to return
% objective values to proper sign
I_max = find(in.maximization == true);
od.fit_save(I_max,:,:) = -od.fit_save(I_max,:,:);
od.fit(I_max,:,:) = -od.fit(I_max,:,:);
od.fit_arch(I_max,:,:) = -od.fit_arch(I_max,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization Stats %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Output termination conditions
if od.iter == in.max_iter
  
  fprintf('\nPareto Solution not necessarily converged.\n');
  fprintf('Maximum number of iterations performed (%g iterations)\n',od.iter);
  
else
  
  fprintf('\n\nPareto Solution Converged\n');
  fprintf('Number of Iterations: %g\n',od.iter);
  
end

% Execution time
fprintf('Time elapsed: %g seconds \n\n',od.time);

% Save results to binary file
save('results.mat','od','in','ui');

return

