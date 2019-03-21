function [] = output(in,od,ui)
%
% [] = output(in,od)
%
% This function is used to allocate memory locations for matrices and
% initialize scalar values.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during 
%                optimization
%
% Output:
%   None. Although no values are returned to a calling function, output is
%   printed to the screen and the results are stored in results.mat.
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitness Conversion %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Determine iteration with best individual - choose first instance
I = find(od.best_fit == min(od.best_fit));
I = min(I);

% Determine individual in iteration with optimal value - choose first instance
p = od.best_member(1,I);

% If performing maximization, objective values were negated. Need to return
% objective values to proper sign
I_max = find(in.maximization == true);
od.fit_save(I_max,:,:) = -od.fit_save(I_max,:,:);
od.fit(I_max,:,:) = -od.fit(I_max,:,:);
od.best_fit(I_max,:) = -od.best_fit(I_max,:);

%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization Stats %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Execution time
fprintf('Time elapsed: %g seconds \n\n',od.time);

%%%%%%%%%%
%% Data %%
%%%%%%%%%%

% Output to screen
fprintf('\nBest Fitness Attained in Iteration %d out of %d',I,od.iter);
fprintf('\nBest Fitness Value: %g',od.best_fit(I));
fprintf('\nTrade Space Location:\n');
fprintf('      %g\n',od.pop(:,p,I));
fprintf('\n\n');

% Save results to binary file
save('results.mat','od','in','ui');

return

