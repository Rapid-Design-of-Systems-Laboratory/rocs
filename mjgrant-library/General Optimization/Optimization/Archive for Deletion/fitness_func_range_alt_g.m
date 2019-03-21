function [fitness,constr] = fitness_func(mat_input,in)
%
% [fitness] = fitness_func(mat_input)
%
% This function is used to return the fitness of all of the test solutions.
% A simulation could be performed in this function or a continuous function
% can be tested. The output returned will represent all of the multi-objective
% fitnesses associated with each test solution. Thus, there is no single-best
% solution.
%
% Input:
%   mat_input - Matrix that assumes that each row specifies a dimension that is 
%               varied to find the optimal solution. Thus, each column 
%               represents the location of the individual.
%
% Output:
%   fitness - matrix containing all of the multi-objective fitness values for 
%             each individual
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Feb. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Reference Trajectories %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete files
!rm -rf dlock
!rm -rf Sim_Info
!mkdir Sim_Info

% Determine number of individuals
num_sims = length(mat_input(1,:));

% Create proper monte_carlo.mc file
write_mci('ref',num_sims-1);

%% Create input data file
fid = fopen('monte_carlo.sts','w');
fprintf(fid,'%12.10e %12.10e %12.10e\n',mat_input);
fclose(fid);

% Create proper monte_tpl.inp file and stress case files
input_05_22('ref_optim');
input_05_22('alt');
input_05_22('steep');
input_05_22('shallow');

% Run reference and stress case 'Monte Carlo'
exec_sims(num_sims,'mc',[],0,'post_msl');

% Determine results
for counter = 1 : 1 : num_sims
  fitness(:,counter) = load(['Sim_Info/b',int2str(counter),'_info.dat']);
end
fitness(1,:) = -fitness(1,:);

% Assign constraint values
constr = fitness(1:2,:);

% Switch fitness values if performing minimization
if in.maximization

  % Change fitness signs
  fitness = -fitness;

end

return

