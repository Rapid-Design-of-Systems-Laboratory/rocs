function run_optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine how function is used %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 0
  clear all; close all; clc;
end

%%%%%%%%%
% Input %
%%%%%%%%%

% Optimization parameters
ub = [20 20]'; % Upper bound
lb = [-20 -20]'; % Lower bound
x0 = [1 1]'; % Initial guess
options = []; % Optimization options

%%%%%%%%%%%%%%%%%%%%%%
%% Run Optimization %%
%%%%%%%%%%%%%%%%%%%%%%

[x,fval] = fmincon('fitness_func_gradient',x0,[],[],[],[],lb,ub,[],options);

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

% Matlab's optimization performance
fprintf('\n\n');
fprintf('Optimization Results:\n');
fprintf('Fitness:  %g\n',-fval);
fprintf('Solution: x=%g, y=%g\n',x(1),x(2));

fprintf('\n');




