function [od] = advance_pop(in,od)
%
% [od] = advance_pop(in,od)
%
% This function advances or initializes the population. The population is
% initialized randomly throughout the entire trade space. The velocity of
% the individuals is also computed by taking into account the cognitive and
% social influences of best solutions.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants for Iteration %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[od] = set_constants(in,od);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update/initialize velocity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[od] = set_velocity(in,od);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update/initialize position %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[od] = set_position(in,od);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = set_constants(in,od)
%
% [od] = set_constants(od,in)
%
% This function assigns the three coefficients used in the velocity
% calculations:
%   r1 - Cognitive parameter
%   r2 - Social parameter
%   w  - Inertia
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Velocity Constant Values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set values for velocity computations that are constant for each iteration
if od.iter == 1
  
  % Initializing
  od.r1(od.iter) = NaN; % Cognitive parameter
  od.r2(od.iter) = NaN; % Social parameter
  od.w(od.iter) = NaN; % Inertia
  
else
  
  % Advancing population
  od.r1(od.iter) = rand; % Cognitive parameter
  od.r2(od.iter) = rand; % Social parameter
  [od] = set_inertia(in,od); % Inertia
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = set_inertia(in,od)
%
% [od] = set_inertia(od,in);
%
% This function calculates the inertia associated with the iteration.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%
%% Compute inertia %%
%%%%%%%%%%%%%%%%%%%%%

% Inertia computed based on type of inertia control
switch in.w_type

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Coefficient of Variance %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Inertia is decremented in a step process when the coefficient of variance is
  % below a given tolerance.
  
  % % Percentage of top crowdedness distances taken for COV calculation
  % in.cov_per = 1.0;
  % 
  % % Coefficient of variance tolerance to modify inertia
  % in.cov_tol = 0.5;
  % 
  % % Inertia geometric factor (<1)
  % in.w_factor = 0.975;
  % 
  % % Initial inertia value
  % in.w_init = 1.4;
  
  case 'cov'

    % Determine calculated crowdedness distance values
    [R,C] = find(isnan(od.cd(1,:,od.iter-1)) == 0);
    C = unique(C);
    
    % Obtain top crowdedness distance values
    [cd_sorted,I] = sort(od.cd(1,C,od.iter-1));
    
    % Remove infinity values for endpoints
    I = find(cd_sorted ~= inf);
    cd_sorted = cd_sorted(I);
    
    % Construct set
    cd_set = cd_sorted(1:ceil(length(I)*in.cov_per));

    % Determine cofficient of variance (COV) of fitness values
    sigma = std(cd_set); % Standard deviation
    mu = mean(cd_set); % Mean
    od.cov(od.iter) = abs(sigma/mu); % COV

    % Initialize inertia, if necessary
    if od.iter == 2
      od.w(od.iter-1) = in.w_init;
    end

    % If coefficient of variance below tolerance, then assume solution 
    % converging ==> Decrement inertia
    if od.cov(od.iter) < in.cov_tol
      od.w(od.iter) = od.w(od.iter-1)*in.w_factor;
    else
      od.w(od.iter) = od.w(od.iter-1);
    end
    
  %%%%%%%%%%%%%%%%%%
  %% Archive Size %%
  %%%%%%%%%%%%%%%%%%
  
  % If archive is full, then assume want to diversify Pareto front. Reduce 
  % inertia to better explore the Pareto front
  
  case 'arch'
    
    % Initialize inertia if necessary
    if od.iter == 2
      od.w(od.iter-1) = in.w_init;
    end
    
    % Determine if archive full
    if ~isnan(od.cd(1,:,od.iter-1))

      % Archive full, decrement inertia
      od.w(od.iter) = od.w(od.iter-1)*in.w_factor;

    else

      % Archive not full, do not modify inertia
      od.w(od.iter) = od.w(od.iter-1);

    end
    
    % Decrement inertia if premature convergence commanded in previous iteration
    if od.flag(od.iter-1) == 1
      od.w(od.iter) = 0;
    end
    
  %%%%%%%%%%%%
  %% Random %%
  %%%%%%%%%%%%

  % Randomly vary inertia between bounds. Usually set around 0.5 for MOPSO
  % search.

  case 'rand'

    % Check to make sure user input proper data
    if in.ib(2) < in.ib(1)

      error('in.ib must be monotonically increasing.');

    end

    % Randomly assign inertia value between upper and lower bounds
    od.w(od.iter) = in.ib(1) + (in.ib(2) - in.ib(1))*rand;
    
  %%%%%%%%%%%%%%%%%%%%%%%%
  %% Migration Distance %%
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  % Compute migration distance of dominant solutions. If maximum distance is
  % less than the threshold input by the user, reduce the inertia.
  
  case 'md'
    
    % Initialize inertia if necessary
    if od.iter == 2
      od.w(od.iter-1) = in.w_init;
    end
    
    % Initialize migration distance
    od.md(od.iter) = 0;
    
    % If enough iterations have passed to compare archives, commence computation
    % of migration distance
    if od.iter > in.arch_del_iter+1
    
      % Determine portion of archives with data
      [R,C] = find(isnan(od.fit_arch(:,:,od.iter-in.arch_del_iter-1)) == 0);
      C_old = unique(C); % Make sure no repetitions
      [R,C] = find(isnan(od.fit_arch(:,:,od.iter-1)) == 0);
      C_new = unique(C); % Make sure no repetitions
      
      % Obtain Pareto front information
      pareto_old = od.fit_arch(:,C_old,od.iter-in.arch_del_iter-1);
      pareto_new = od.fit_arch(:,C_new,od.iter-1);

      [row_old,col_old] = size(pareto_old);

      % Determine which points in new Pareto front dominate points in old Pareto
      % front. Then, calculate the migration distance the dominant points have
      % moved.

      % Determine portion of archive with data
      C = C_new;

      % Loop through each individual to test dominance
      for ctr = 1 : 1 : length(C)

        % Create matrix repeating fitness values for individual (for matrix math)
        index = C(ctr);
        fit = od.fit_arch(:,index,od.iter-1);
        fit_new = fit(:,ones(size(pareto_old,2),1));

        % Create matrix of relative fitness values
        diff_set = pareto_old - fit_new;

        % Determine characteristics of column in augmented matrix (for dominance)
        [r_big,c_big] = find(diff_set > 0);
        [r_small,c_small] = find(diff_set < 0);
        c_small = unique(c_small);
        c_big = unique(c_big);

        % Determine dominance
        if length(c_big) == col_old

          % Individual is not dominated by any other individual in old Pareto
          % front. If individual is in both archives, then is part of
          % nondominated solutions and if statement will not execute.
          
          % Determine which members of the old Pareto front are dominated
          dominates = setdiff(c_big,c_small)';

          if ~isempty(dominates)
            
            % Compute migration distance of dominant member in archive
            rel_pos = fit_new(:,dominates) - pareto_old(:,dominates);
            
            for ctr_pos = 1 : 1 : length(dominates)
              
              dist = norm(rel_pos(:,ctr_pos));
              if dist > od.md(od.iter)
                od.md(od.iter) = dist;
              end
              
            end
            
          end
          
        end

      end
      
      % If no new dominant members, migration distance is set to zero
      if isnan(od.md(od.iter))
        
        od.md(od.iter) = 0;
        
      end

      % If normalized migration distance less than user input, reduce inertia
      if od.md(od.iter) <= in.dist
        od.w(od.iter) = od.w(od.iter-1)*in.w_factor;
      else
        od.w(od.iter) = od.w(od.iter-1);
      end
      
    else
      
      % Not enough iterations have passed
      od.w(od.iter) = od.w(od.iter-1);
      
    end

  otherwise

    error('Invalid in.w_type value.');
    
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = set_velocity(in,od)
%
% [od] = set_velocity(in,od)
%
% This function computes the velocity for each individual. This determines how
% the population moves about the trade space.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%
%% Set Velocity %%
%%%%%%%%%%%%%%%%%%

% Determine if first iteration
if od.iter == 1
  
  % Randomly initialize velocity
  rand_mat = rand(od.dim,in.num_ind);
  diff_mat = in.vel_max - in.vel_min;
  od.vel(:,:,od.iter) = rand_mat.*diff_mat(:,ones(size(rand_mat,2),1)) + ...
    in.vel_min(:,ones(size(rand_mat,2),1));

else
  
  % Determine archive member for social influence of Pareto front. Determine
  % portion of archive with least crowded members.
  % Determine crowdedness distance order in archive
  [Y,I_crowd] = sort(od.cd(1,:,od.iter-1),2,'descend');
  [R,C] = find(isnan(Y) == 0);
  C = unique(C);
  I_crowd = I_crowd(C);
  
  % Determine upper bound for social influence selection from archive
  I_upper = floor(length(C)*in.cp);
  
  % Ensure upper bound for social influence selection from archive includes the
  % two boundary points
  if I_upper < 2
    
    % At least include two boundary points that have Inf crowdedness distance
    I_upper = 2;
    
  end
  
  % If have one non-dominated solution, objectives may not be in competition.
  % Have the swarm follow the leader
  if length(C) == 1
    
    fprintf('  Warning! Only one non-dominated solution!');
    I_upper = 1;
    
  end
  
  % Cycle through each individual to update velocity. Randomly select social 
  % term from archive based on having high crowdedness distance.
  % Randomly select individuals from social portion of archive
  I_select = ceil(rand(1,in.num_ind)*I_upper);

  % Determine location of individuals in history matrix
  n = I_crowd(I_select);

  %%%%%%%%%%%%%%%%%%%%%
  %% Update Velocity %%
  %%%%%%%%%%%%%%%%%%%%%

  % Determine each term in the velocity equation
  term1 = od.w(od.iter)*od.vel(:,:,od.iter-1);
  term2 = in.c1*od.r1(od.iter)*(od.best(:,:,od.iter-1) - od.pop(:,:,od.iter-1));
  term3 = in.c2*od.r2(od.iter)*(od.pos_arch(:,n,od.iter-1) - ...
    od.pop(:,:,od.iter-1));

  % Update velocity
  od.vel(:,:,od.iter) = in.X*(term1 + term2 + term3);

  %%%%%%%%%%%%%%%%%%%%%%%%
  %% Constrain Velocity %%
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  % Determine size of velocity matrix for indexing bias
  [row,col] = size(od.vel(:,:,od.iter));
  
  % Determine indexing bias based on third dimension (iteration number)
  index_bias = row*col*(od.iter-1);
  
  % Construct matrix of min vel that is the same size as the vel matrix
  temp_vel = in.vel_min(:,ones(size(od.vel,2),1));
  
  % Determine indices of velocity components that exceed constraints
  I = find(od.vel(:,:,od.iter) - temp_vel < 0);
  
  % Force velocity components to be at boundary of constraint
  od.vel(I+index_bias) = temp_vel(I);
  
  % Construct matrix of max vel that is the same size as the vel matrix
  temp_vel = in.vel_max(:,ones(size(od.vel,2),1));
  
  % Determine indices of velocity components that exceed constraints
  I = find(od.vel(:,:,od.iter) - temp_vel > 0);
  
  % Force velocity components to be at boundary of constraint
  od.vel(I+index_bias) = temp_vel(I);

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = set_position(in,od)
%
% [od] = set_position(in,od)
%
% This function advances/initializes the population based on the velocity
% calculations.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during optimization
%
% Output:
%   od - structure containing variables of interest during optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%
%% Set Position %%
%%%%%%%%%%%%%%%%%%

% Determine if first iteration
if od.iter == 1
  
  % Randomly initialize position (potentially covering entire trade space)
  rand_mat = rand(od.dim,in.num_ind);
  diff_mat = in.pos_max - in.pos_min;
  od.pop(:,:,od.iter) = rand_mat.*diff_mat(:,ones(size(rand_mat,2),1)) + ...
    in.pos_min(:,ones(size(rand_mat,2),1));
  
  % Determine infeasible solutions and re-initialize infeasible solutions until
  % have entire population with a feasible initial guess
  
  if od.num_constr > 0
    
    % Compute fitness of population
    [od.fit_save(:,:,od.iter),od.constr(:,:,od.iter)] = eval_pop(od.pop(:,:,od.iter));

    % Determine infeasible solutions and modify od.fit to NaNs
    I_broke = [];
    for ctr = 1 : 1 : od.num_constr
      I_satisfy = eval(['find(od.constr(ctr,:,od.iter)',in.constr{ctr,2},num2str(in.constr{ctr,3}),')']);
      I_broke = [I_broke setdiff([1:1:in.num_ind],I_satisfy)];
    end
    C = unique(I_broke);
    od.fit(:,:,od.iter) = od.fit_save(:,:,od.iter);
    od.fit(:,C,od.iter) = NaN;

    % Change infeasible solutions to feasible (randomly)
    % Iterate until have all feasible solutions
    while ~isempty(C)

      % Output to screen
      fprintf(['%g infeasible solutions. Reinitialize infeasible ', ...
        'solutions.\n'],length(C));

      % Modify initial guess of infeasible solutions
      rand_mat = rand(od.dim,length(C));
      diff_mat = in.pos_max - in.pos_min;
      od.pop(:,C,od.iter) = rand_mat.*diff_mat(:,ones(size(rand_mat,2),1)) + ...
        in.pos_min(:,ones(size(rand_mat,2),1));

      % Evaluate fitness of new solutions
      [od.fit_save(:,C,od.iter),od.constr(:,C,od.iter)] = eval_pop(od.pop(:,C,od.iter));

      % Determine infeasible solutions and modify od.fit to NaNs
      I_broke = [];
      for ctr = 1 : 1 : od.num_constr
        I_satisfy = eval(['find(od.constr(ctr,:,od.iter)',in.constr{ctr,2},num2str(in.constr{ctr,3}),')']);
        I_broke = [I_broke setdiff([1:1:in.num_ind],I_satisfy)];
      end
      C = unique(I_broke);
      od.fit(:,:,od.iter) = od.fit_save(:,:,od.iter);
      od.fit(:,C,od.iter) = NaN;

    end
    
  else
    
    % Compute fitness of population
    [od.fit_save(:,:,od.iter)] = eval_pop(od.pop(:,:,od.iter));
    od.fit(:,:,od.iter) = od.fit_save(:,:,od.iter);
    
  end
  
  % Turn maximization objectives into minimization
  I = find(in.maximization == true);
  od.fit_save(I,:,od.iter) = -od.fit_save(I,:,od.iter);
  od.fit(I,:,od.iter) = -od.fit(I,:,od.iter);
    
else
  
  %%%%%%%%%%%%%%%%%%%%%
  %% Update Position %%
  %%%%%%%%%%%%%%%%%%%%%
  
  % Move each individual in population
  od.pop(:,:,od.iter) = od.pop(:,:,od.iter-1) + od.vel(:,:,od.iter);
  
  % Determine size of velocity matrix for indexing bias
  [row,col] = size(od.pop(:,:,od.iter));
  
  % Determine indexing bias based on third dimension (iteration number)
  index_bias = row*col*(od.iter-1);

  % Determine if upper bound exceeded. If so, place particle at wall
  temp_pos = in.pos_max(:,ones(size(od.pop,2),1));
  I_upper = find(od.pop(:,:,od.iter) > temp_pos);
  od.pop(I_upper+index_bias) = temp_pos(I_upper);

  % Determine if lower bound exceeded. If so, place particle at wall
  temp_pos = in.pos_min(:,ones(size(od.pop,2),1));
  I_lower = find(od.pop(:,:,od.iter) < temp_pos);
  od.pop(I_lower+index_bias) = temp_pos(I_lower);

  % Listing of changed dimensions
  I = union(I_upper,I_lower);

  % Bounce or set vel to zero for vel components that exceed constraint
  value = floor(2*rand(1,length(I))); % Randomly choose bounce/rest
  
  % Bounce
  I_bounce = find(value >= 1);
  od.vel(I_bounce+index_bias) = -od.vel(I_bounce+index_bias);
  
  % Set velocity to zero
  I_stop = find(value < 1);
  od.vel(I_stop+index_bias) = 0;

end

return
