function [od] = post_fitness_calcs(in,od)
%
% [od] = post_fitness_calcs(in,od)
%
% This function is used to perform any calculations associated with the
% current iteration that require knowledge of the fitnesses of the
% individuals. This function is executed before the termination criteria is
% checked.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during 
%                optimization
%
% Output:
%   od - structure containing variables of interest during 
%                optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


% Determine if particles moved to a higher fitness solution in previous 
% iteration
[od] = better_fitness_calc(od);

% Determine location of best neighbors
[od] = neighbor_calcs(in,od);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = neighbor_calcs(in,od)
%
% [od] = neighbor_calcs(od,in)
%
% This function is used to determine the individuals that serve as the most fit
% neighbors.
%
% Input:
%   in - structure containing inputs from input deck
%   od - structure containing variables of interest during 
%                optimization
%
% Output:
%   od - structure containing variables of interest during 
%                optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Location of Best Neighbor %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through each individual in population
for p = 1 : 1 : in.num_ind
  
  % Set current particle position vector
  pos = od.pop(:,p,od.iter);
  
  % Duplicate vector into matrix for matrix subraction
  pos = pos(:,ones(size(od.pop(:,:,od.iter),2),1));
  
  % Determine distance each individual is from test individual
  dist = sqrt(sum((pos - od.pop(:,:,od.iter)).^2,1));
  
  % Sort distances
  [temp,I] = sort(dist);
  
  % Determine neighbors
  neighbors = I(2:in.num_nbr+1)';
  
  % Determine best neighbor of particle
  fit = od.bestSoFar(1,neighbors,od.iter);
  [val,I] = min(fit);
  
  % Best neighbor location
  % Determine if a feasible neighbor exists
  if isnan(val)
    
    % No feasible neighbor, copy best neighbor from previous iteration
    od.best_nbr(:,p,od.iter) = od.best_nbr(:,p,od.iter-1);
    
  else
    
    % Feasible neighbor exists
    od.best_nbr(:,p,od.iter) = od.pop(:,neighbors(I),od.iter);    
    
  end
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = better_fitness_calc(od)
%
% [od] = better_fitness_calc(od)
%
% This function is used to determine which individual moved to a solution that 
% is better than the individual has seen during the entire optimization process.
%
% Input:
%   od - structure containing variables of interest during 
%                optimization
%
% Output:
%   od - structure containing variables of interest during 
%                optimization
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Particles Best Memory %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if first iteration
if od.iter ~= 1
  
  % Determine which individuals have found new best solution
  od.bestSoFar(1,:,od.iter) = od.bestSoFar(1,:,od.iter-1);
  od.best(:,:,od.iter) = od.best(:,:,od.iter-1);
  I = find(od.bestSoFar(1,:,od.iter-1) - od.fit(1,:,od.iter) > 0);
  
else
  
  % First iteration, every individual is located at the best solution
  I = 1 : 1 : length(od.fit(1,:,od.iter));
  
end

% Assign data to matrices recording the best solutions
od.bestSoFar(1,I,od.iter) = od.fit(1,I,od.iter);
od.best(:,I,od.iter) = od.pop(:,I,od.iter);
od.best_fit(1,od.iter) = min(od.fit(1,:,od.iter));

return
