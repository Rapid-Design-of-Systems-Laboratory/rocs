function [] = ea_optim()
%
% [] = ea_optim()
%
% This function executes the proper evolutionary optimization algorithm 
% specified by the user in set_input.m. Functions corresponding to the
% specified algorithm are located in the correpsonding directory. The user must
% supply a function to evaluate the objectives and any constraints (located in
% eval_obj.m (or the default eval_obj.m in the directory corresponding to the 
% optimization algorithm is used). Output is returned in results.mat located
% in the working directory. Throughout the optimization, a restart.mat file may
% be created if specified by the user. This file is identical to results.mat but
% is a snapshot of the current optimization.
%
% Input:
%   None. Inputs specified in set_input.m
%
% Output:
%   None. Results saved in results.mat (and restart.mat)
%
% Optimization algorithms supported:
%   Particle Swarm Optimization (located in the 'pso' directory)
%   Genetic Algorithm (located in the 'ga' directory)
%   Multi-objective Particle Swarm Optimization (located in the 'mopso'
%     directory
%
% Note: Additional algorithms must follow the template in ea_optim.m. Note that
% the function names are consistent for each optimization algorithm but reside 
% in the directory corresponding to the optimization algorithm.
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain inputs/initialize %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear everything if not compiled. Otherwise would get message saying clc is
% not required for compiled executables.
if ~isdeployed
  clc;
end

% Get time and date of beginning of optimization
t0 = clock;

% Determine if restarting (look for restart file)
if exist('./restart.mat','file') % Restart
  
  % Load data
  temp = load('./restart.mat');
  in = temp.in;
  od = temp.od;
  ui = temp.ui;
  
  % Tell user using restart.mat file
  fprintf('Using restart.mat file (from iter %g)\n',od.iter);
  
else % Starting from new

  % Set required variables before obtaining input
  od.iter = 1; % Iteration number
  
  % Save eval_obj.m information
  fid = fopen('eval_pop.m','r');

  % Read eval_obj.m file
  index = 0;
  line = fgets(fid);
  while line ~= -1
    index = index + 1;
    ui.eval_pop{index} = line;
    line = fgets(fid); % Read next line
  end

  % Close eval_obj.m
  fclose(fid);
  
  % Save set_input.m information
  fid = fopen('set_input.m','r');

  % Read set_input.m file
  index = 0;
  line = fgets(fid);
  while line ~= -1
    index = index + 1;
    ui.set_input{index} = line;
    line = fgets(fid); % Read next line
  end

  % Close set_input.m
  fclose(fid);

  % Obtain input from set_input.m
  [in] = set_input;

  % Initialize variables and matrices - algorithm specific
  [od] = initialize(in,od);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop of optimizer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while true
  
  % Advance/initialize population
  [od] = advance_pop(in,od);  
  
  % Compute objectives and constraints of population. For algorithms that 
  % support constrained optimization, the objectives may already be computed for
  % the first iteration. If this is the case, do not waste time re-evaluating 
  % the objectives. ea_optim.m assumes that the objective values are initialized
  % as NaNs.
  if isnan(od.fit(:,:,od.iter))
    
    % Objectives not already computed, evaluate objectives (and constraints if
    % applicable).
    if od.num_constr > 0
      
      % Constrained problem. Query eval_pop.m with 2 outputs (objective
      % information and constraint information).
      [od.fit_save(:,:,od.iter),od.constr(:,:,od.iter)] = eval_pop(od.pop(:,:,od.iter));
      % Determine individuals that broke constraint and set to NaN
      I_broke = [];
      for ctr = 1 : 1 : od.num_constr
        I_satisfy = eval(['find(od.constr(ctr,:,od.iter)',in.constr{ctr,2},num2str(in.constr{ctr,3}),')']);
        I_broke = [I_broke setdiff([1:1:in.num_ind],I_satisfy)];
      end
      I_broke = unique(I_broke);
      % Assign NaNs to individuals who broke constraints
      od.fit(:,:,od.iter) = od.fit_save(:,:,od.iter);
      od.fit(:,I_broke,od.iter) = NaN;
      
    else
      
      % Unconstrained problem. Query eval_pop.m with only 1 output (objective
      % information).
      [od.fit_save(:,:,od.iter)] = eval_pop(od.pop(:,:,od.iter));
      od.fit(:,:,od.iter) = od.fit_save(:,:,od.iter);
      
    end
    
    % Turn maximization objectives into minimization. All algorithms are assumed
    % to be minimizing the objectives. This is only used internally. The results
    % are returned to their normal values in output.m.
    I = find(in.maximization == true);
    od.fit_save(I,:,od.iter) = -od.fit_save(I,:,od.iter);
    od.fit(I,:,od.iter) = -od.fit(I,:,od.iter);
    
  end
  
  % Perform calculations after objectives and constraints have been evaluated
  % (if necessary)
  [od] = post_obj_calcs(in,od);

  % Real-time optimization information output to screen to allow user to monitor
  % progress of the optimization.
  fprintf('Iter: %g/%g', od.iter,in.max_iter);
  if strcmp(in.method,'pso') || strcmp(in.method,'mopso')
    fprintf(' / w: %g',od.w(od.iter));
    if strcmp(in.method,'pso')
      fprintf(' / BF: %g / COV: %g',od.best_fit(od.iter),od.cov(od.iter));
    elseif strcmp(in.method,'mopso')
      % Determine size of archive
      [R,C] = find(isnan(od.fit_arch(:,:,od.iter)) == 0);
      C = unique(C);
      fprintf(' / Arch: %g/%g',length(C),in.num_arch);
    end
  end
  fprintf('\n');
  
  % Termination criteria for optimization. Information from this function may be
  % passed in od structure to influence other functions of the optimizer (namely
  % flags denoted premature termination was commanded).
  [od,terminate] = check_terminate(in,od);
  if terminate
    break;
  end
  
  % Perform intermediate saving (save optimization information to restart.mat)
  if mod(od.iter,in.freq_save) == 0
    save('restart.mat','od','in','ui');
  end

  % Step iteration
  od.iter = od.iter + 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Final calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Compute execution time. This is performed using the difference in date and
% time. This information could be inaccurate if a restart file was used (thus
% adding time in which the algorithm was not running).
t1 = clock;
od.time = etime(t1,t0);

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%
  
% Execute output.m
output(in,od,ui);

return

