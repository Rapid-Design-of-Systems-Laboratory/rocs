function [od] = advance_pop(in,od)
%
% [od] = advance_pop(in,od)
%
% This function advances or initializes the population. The population is
% initialized randomly throughout the entire trade space.
%
% Population advancing processess:
%   Parent selection - Random process but fit individuals have a higher chance 
%                      of being selected
%   Children formation - Genes from parents blended via crossover
%   Children mutation - Genes randomly changed to potentially test solutions in 
%                       new regions of the trade space
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during optimization
%
% Output:
%   od - structure containing variables of interest during optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advance/Initialize Population %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform check if population must be initialized
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
  % Determine parents %
  %%%%%%%%%%%%%%%%%%%%%
  
  [od] = determine_parents(in,od);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Determine children - Perform crossover %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [od] = determine_children(in,od);
  
  %%%%%%%%%%%%%%%%%%%
  % Mutate children %
  %%%%%%%%%%%%%%%%%%%
  
  [od] = mutate_children(in,od);
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = determine_parents(in,od)
%
% [od] = determine_parents(in,od)
%
% This function forms the mating pool by randomly selecting individuals to mate
% using a roulette wheel approach. Individuals with fitness values have an
% increased chance to mate. If elitism is on, the best individual is directly 
% copied to the next generation and must be assigned as a parent.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during optimization
%
% Output:
%   od - structure containing variables of interest during optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Make all fitness values positive and start at zero for roulette selection
fitness_abs = -(od.fit(1,:,od.iter-1) - max(od.fit(1,:,od.iter-1)));
sum_fitness_abs = sum(fitness_abs);
cum_fitness_abs = cumsum(fitness_abs);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spin Roulette Wheel %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine pointers for roulette wheel spinning. Need to randomly select
% parents to equal original population size. Elite member already selected as 
% parent if elitism on.
pointer = sum_fitness_abs*rand(1,in.num_ind);

% Determine if elitism on for parent selection
if in.elitism
  
  % Add elite member as parent
  od.parents(:,1,od.iter) = od.pop(:,od.best_member(od.iter-1),od.iter-1);

end

% Determine which members of population the pointers are identifying and add
% individual to parent selection pool. Skip over first pointer if elite member
% already selected as parent.
for counter = 1+in.elitism : 1 : in.num_ind
  
  % Determine pointer location
  I = find(pointer(counter) <= cum_fitness_abs);
  I = min(I);

  % Add parent to parent selection pool
  od.parents(:,counter,od.iter) = od.pop(:,I,od.iter-1);
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = determine_children(in,od)
%
% [od] = determine_children(in,od)
%
% This function performs reproduction of the parents using crossover. The point 
% for crossover is randomly chosen for each mating pair.
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
%% Elite Member %%
%%%%%%%%%%%%%%%%%%

% Transfer elite member if necessary
if in.elitism
  
  % Copy elite member into next geneation (elite member is first parent)
  od.pop(:,1,od.iter) = od.parents(:,1,od.iter);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Form Children - Crossover %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p1 is individual who gets to mate
for p1 = 1+in.elitism : 1 : in.num_ind

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Determine Mate Randomly %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Find a random mate other than itself
  p2_set = setdiff(1:1:in.num_ind,p1); % Second parent set (p1 removed)
  p2 = p2_set(ceil(rand*length(p2_set))); % Select second parent
  
  %%%%%%%%%%%%%%%%%%%%%%%
  %% Perform Crossover %%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  % Determine if crossover should be performed
  if in.cross_prob > rand
    
    % Perform crossover by randomly choosing site for crossover
    site = ceil(rand*od.dim);
    
    % Ensure crossover point is not last trait
    if site == od.dim
      site = site - 1;
    end
    
    % Ensure site is not last trait (otherwise have direct copy)
    if site == od.dim
      
      % Reduce site to second to last trait
      site = od.dim - 1;
      
    end
    
    % Form the child by the swapping of genetic material between the parents
    od.pop(:,p1,od.iter) = [od.parents(1:site,p1,od.iter); ...
      od.parents(site+1:od.dim,p2,od.iter)];
    
  else
    
    % No crossover occurs
    % Copy non-crossovered chromosomes into next generation
    % In this case, simply take first parent and make it the child.
    od.pop(:,p1,od.iter) = od.parents(:,p1,od.iter);
    
  end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = mutate_children(in,od)
%
% [od] = determine_children(in,od)
%
% This function randomly mutates children in order to enhance global searching
% capability.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during optimization
%
% Output:
%   od - structure containing variables of interest during optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%
%% Mutate Children %%
%%%%%%%%%%%%%%%%%%%%%

% Loop through each child. Do not mutate elite individual. Elite member is first
% member in population (if elitism on).
for pop_member = 1+in.elitism : 1 : in.num_ind
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Determine if Individual Mutates %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Probability of individual mutating based on input
  if in.mutate_prob_ind > rand

    %%%%%%%%%%%%%%%%%%%
    %% Mutate Traits %%
    %%%%%%%%%%%%%%%%%%%
    
    % Determine which traits to mutate
    mutate_traits = find(in.mutate_prob_trait > rand(1,od.dim));
    
    % Create random traits - Technically could be same as before
    rand_traits = rand(length(mutate_traits),1).*(in.pos_max(mutate_traits) ...
      - in.pos_min(mutate_traits)) + in.pos_min(mutate_traits);
    
    % Make sure traits different than before. If so, alter all traits
    while length(setdiff(od.pop(mutate_traits,pop_member,od.iter),...
        rand_traits)) ~= length(rand_traits)
      rand_traits = rand(length(mutate_traits),1).* ...
        (in.pos_max(mutate_traits) - in.pos_min(mutate_traits)) + ...
        in.pos_min(mutate_traits);
    end
    
    % Assigned mutated traits to population
    od.pop(mutate_traits,pop_member,od.iter) = rand_traits;
    
  end

end

return

