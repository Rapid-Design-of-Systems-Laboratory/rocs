function [od] = advance_pop(in,od)
%
% [od] = advance_pop(in,od)
%
% This function advances or initializes the population throughout the trade
% space. The population is
% initialized randomly throughout the entire trade space. The velocity of
% the individuals is also computed by taking into account the cognitive and
% social influences of best solutions.
%
% Population initialization process:
%		Position initialization - The position of the population is randomly 
% 	initialized throughout the trade space (bounded by the input of the user).
%		This is used to promote diversity of the initial population.
%
% Population advancing processess:
%   Set velocity - Influenced by the inertia, location of particle's best
%   	solution (cognitive influence), and location of best neighbor (social
%			influence)
%   Move individuals - Perform step in velocity direction for each individual. 
%			If individual would exceed limits of trade space, then the individual is 
%			placed on the boundary and the velocity is modified.
%				Velocity modification - Two different operations are used to modify the
%					velocity of individuals placed at the boundary with both operations
%					having an equal probability of execution. First, the velocity of the
%					individual may be negated (resulting in an effective "bouncing" off
%					the boundary). Second, the individual might be placed at rest. Both
%					methods are utilized due to the varying effectiveness of each method
%					for various problems.
%
% Input:
%   in - structure containing inputs from input deck (set_input.m)
%   od - structure containing variables of interest during optimization
%
% Output:
%   od - structure containing variables of interest during optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Values Constant for Iteration %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% calculations to the od structure:
%   r1 - cognitive parameter
%   r2 - social parameter
%   w  - inertia
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during optimization
%
% Output:
%   od - structure containing variables of interest during optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Velocity Constants %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if first iteration for initialization
if od.iter == 1
  % Initializing. Set constants to NaN.
  od.r1(od.iter) = NaN; % Cognitive parameter
  od.r2(od.iter) = NaN; % Social parameter
  od.w(od.iter) = NaN; % Inertia
else
  % Advancing population
  od.r1(od.iter) = rand; % Cognitive parameter (random value - (0,1))
  od.r2(od.iter) = rand; % Social parameter (random value - (0,1))
  [od] = set_inertia(in,od); % Inertia
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = set_inertia(in,od)
%
% [od] = set_inertia(in,od);
%
% This function calculates the inertia associated with the iteration. Initially,
% the inertia is set "high" to promote global searching and is reduced as the
% iterations progress to promote local searching. An inertia above 1 corresponds
%	to acceleration of particles throughout the trade space. An inertia below 1
%	corresponds to a velocity decay. A value of 1.4 is considered quite high for
%	an initial global search and should transition near zero for local searching. 
%	In certain cases, the inertia is set much above 1 since several iterations may
%	be performed before the statistics build to a sufficient level. During these
%	iterations, the inertia may be geometrically reduced several times, and the 
%	algorithm should remain in a global searching mode when the statistics build
%	to the sufficient level.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during optimization
%
% Output:
%   od - structure containing variables of interest during optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%
%% Compute Inertia  %%
%%%%%%%%%%%%%%%%%%%%%%

% Inertia computed based on type of decay
switch in.w_type

  %%%%%%%%%%%%
  %% Linear %%
  %%%%%%%%%%%%
  
  % Inertia linearly decays to a value = initial inertia - 1
  
  case 'linear'
    
    od.w(od.iter) = in.w_init - od.iter/in.max_iter;

  %%%%%%%%%%%%%%%
  %% Quadratic %%
  %%%%%%%%%%%%%%%
  
  % Inertia quadratically decays from initial value to zero at max iterations
  
  case 'quadratic'

    od.w(od.iter) = in.w_init/(1-in.max_iter)^2*(od.iter-in.max_iter)^2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Coefficient of Variance %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Inertia is decremented in a step process when the coefficient of variance is
  % below a given tolerance (denoting the algorithm has sufficiently gained an
  %	understanding of the trade space and increment towards the local search).
  
  case 'cov'

    % Current best fitness set
    fitness = od.bestSoFar(1,:,od.iter-1);

    % Obtain top fitness values
    [fitness_sorted,I] = sort(fitness);
    upper_value = ceil(length(fitness)*in.cov_per);
    fitness_set = fitness_sorted(1:upper_value);

    % Determine cofficient of variance (COV) of top fitness values
    sigma = std(fitness_set); % Standard deviation
    mu = mean(fitness_set); % Mean
    od.cov(od.iter) = abs(sigma/mu); % COV

    % Initialize inertia if necessary
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
% the population moves about the trade space. The velocity of each individual is
%	influenced by the cognitive parameter (the best location in memory), the
%	social parameter (communication with other members of the swarm), and the
%	inertia term (the desire for the individual to continue exploring in the 
%	direction it is traveling).
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
%% Set Velocity %%
%%%%%%%%%%%%%%%%%%

% Determine if first iteration for initialization
if od.iter == 1

  % Randomly initialize velocity within the bounds specified by the user
  rand_mat = rand(od.dim,in.num_ind);
  diff_mat = in.vel_max - in.vel_min;
  od.vel(:,:,od.iter) = rand_mat.*diff_mat(:,ones(size(rand_mat,2),1)) + ...
    in.vel_min(:,ones(size(rand_mat,2),1));

else
  
  %%%%%%%%%%%%%%%%%%%%%
  %% Update Velocity %%
  %%%%%%%%%%%%%%%%%%%%%
  
  % Inertia component of velocity
  term1 = od.w(od.iter)*od.vel(:,:,od.iter-1);
  
  % Cognative (best remembered) component of velocity
  term2 = in.c1*od.r1(od.iter)*(od.best(:,:,od.iter-1) - od.pop(:,:,od.iter-1));
  
  % Social (best neighbor) component of velocity
  term3 = in.c2*od.r2(od.iter)*(od.best_nbr(:,:,od.iter-1) - ...
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
  
  % Construct matrix of min velocity that is the same size as the vel matrix
  temp_vel = in.vel_min(:,ones(size(od.vel,2),1));
  
  % Determine indices of velocity components that exceed min velocity
  I = find(od.vel(:,:,od.iter) - temp_vel < 0);
  
  % Force velocity components to be at boundary of velocity constraints
  od.vel(I+index_bias) = temp_vel(I);
  
  % Construct matrix of max velocity that is the same size as the vel matrix
  temp_vel = in.vel_max(:,ones(size(od.vel,2),1));
  
  % Determine indices of velocity components that exceed max velocity
  I = find(od.vel(:,:,od.iter) - temp_vel > 0);
  
  % Force velocity components to be at boundary of velocity constraint
  od.vel(I+index_bias) = temp_vel(I);
  
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = set_position(in,od)
%
% [od] = set_position(in,od)
%
% This function advances/initializes the population (within the trade space) 
%	based on the velocity calculations.
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

    % Continue to randomly initialize infeasible solutions within the trade space
    %	until feasible.
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

  % Listing of changed dimensions (to enforce trade space constraints)
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

