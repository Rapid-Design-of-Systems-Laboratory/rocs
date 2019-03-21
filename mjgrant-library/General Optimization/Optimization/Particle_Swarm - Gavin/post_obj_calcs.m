function [od] = post_obj_calcs(in,od)
%
% [od] = post_obj_calcs(in,od)
%
% This function is used to perform any calculations associated with the
% current iteration that require knowledge of the objectives and constraints for
% the population. This function is executed before the termination criteria is
% checked.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Feb. 2006


%%%%%%%%%%%%%%%%%%%%
%% Update Archive %%
%%%%%%%%%%%%%%%%%%%%

[od] = modify_archive(in,od);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine New Pareto Optimal Locations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if particles moved to a new non-dominated solution
% Maybe change such that new solution must dominate old one
[od] = higher_fitness_calc(in,od);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = modify_archive(in,od)
%
% [od] = modify_archive(in,od)
%
% This function is used to update the external archive containing the
% Pareto-optimal solutions
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Feb. 2006


%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Copy old archive to next iteration if not first iteration
if od.iter ~= 1
  
  % Initialize fitness archive
  od.fit_arch(:,:,od.iter) = od.fit_arch(:,:,od.iter-1);
  
  % Initialize position archive
  od.pos_arch(:,:,od.iter) = od.pos_arch(:,:,od.iter-1);
  
end

% Determine portion of archive with data
[R,C] = find(isnan(od.fit_arch(:,:,od.iter)) == 0);
C = unique(C); % Make sure no repetitions

% Determine portion of current population that did not violate any constraints.
% Violated individuals are represented as NaNs in objective matrix.
[R,C_fit] = find(isnan(od.fit(:,:,od.iter)) == 0);
C_fit = unique(C_fit); % Make sure no repetitions

% Create augmented matrix combining objectives of individuals (that did not
% violate any constraints) and archive matrix
augmented_matrix = [od.fit(:,C_fit,od.iter) od.fit_arch(:,C,od.iter)];

% Determine size of matrices. col_fit is used to loop through each individual to
% test for Pareto optimality. col_augmented is used for Pareto dominance test
% below.
[row_fit,col_fit] = size(od.fit(:,C_fit,od.iter));
[row_augmented,col_augmented] = size(augmented_matrix);

% Vector of individuals that are non-dominated
non_dominated_set = [];

% Vector of members of the archive that are dominated
dominated_archive_set = [];

% Vector of members in augmented matrix that are dominated by individuals
dominates = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test for Pareto Optimality %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through each individual to test dominance
for fitness_counter = 1 : 1 : col_fit

  % Create matrix repeating objective values for individual (for matrix math)
  fit = od.fit(:,C_fit(fitness_counter),od.iter);
  fit_new = fit(:,ones(size(augmented_matrix,2),1));

  % Create matrix of relative objective values. If entire column is positive,
  % then the individual associated with that column is dominated by individual.
  % If at least one element of a column is positive, then the individuals are
  % non-dominated with respect to each other. If entire column is negative, then
  % individual associated with that column dominates the test individual.
  diff_set = augmented_matrix - fit_new;

  % Determine characteristics of column in augmented matrix (for dominance)
  [r_big,c_big] = find(diff_set > 0);
  [r_small,c_small] = find(diff_set < 0);
  c_small = unique(c_small);
  c_big = unique(c_big);

  % Determine dominance. If at least one element in each column is positive
  % (except for column representing test individual), then individual is
  % non-dominated.
  if length(c_big) == (col_augmented - 1)
    
    % Individual is not dominated by another other individual. Determine which 
    % members of the augmented matrix are dominated. c_small contains columns
    % representing individuals that are better in at least one objective when
    % compared to the test individual. Thus, those individuals that are worse in
    % at least one objective (contained in c_big) but not better in any other
    % objectives (as would be contained in c_small) are dominated. Dominates 
    % includes values from entire augmented matrix.
    dominates = [dominates setdiff(c_big,c_small)'];
    
    % Store individual as non-dominated solution
    non_dominated_set = [non_dominated_set C_fit(fitness_counter)];
    
  end

end

% Dominates includes values from entire augmented matrix. Trim to determine
% those solutions from archive that are dominated.
dominates = dominates - col_fit;
I = find(dominates > 0);
dominated_archive_set = C(unique(dominates(I)));

% "Remove" dominated elements from objective archive (change to NaNs)
od.fit_arch(:,dominated_archive_set,od.iter) = NaN;

% "Remove" elements from position archive (change to NaNs)
od.pos_arch(:,dominated_archive_set,od.iter) = NaN;

% Determine vacant spots in archive
[R_NaN,C_NaN] = find(isnan(od.fit_arch(:,:,od.iter)) == 1);
C_NaN = unique(C_NaN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Elements to Archive %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add elements to archive. Trim archive to proper size if necessary by removing
% elements with lowest crowdedness distance.
for counter = 1 : 1 : length(non_dominated_set)
  
  if counter <= length(C_NaN)
    
    % Archive not full. Add element to open slot.
    C = C_NaN(counter); % Index of element to add
    
  else
    
    % Archive full. Determine which element to replace. Member of archive with
    % snallest crowdedness distance is removed and is replaced by new
    % individual.
    [od] = determine_crowdedness(in,od); % Determine crowdedness distance
    
    % Determine member with smallest crowdedness distance
    [Y,I] = sort(od.cd(1,:,od.iter));
    
    % Assign element to replace
    C = I(1);
    
  end
  
  % Add non-dominated element to archive
  od.fit_arch(:,C,od.iter) = od.fit(:,non_dominated_set(counter),od.iter);

  % Add non-dominated element to position matrix associated with archive
  od.pos_arch(:,C,od.iter) = od.pop(:,non_dominated_set(counter),od.iter);
  
end

% Update crowdedness distance since all new individuals have been added to the
% archive.
[od] = determine_crowdedness(in,od);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = higher_fitness_calc(in,od)
%
% [od] = higher_fitness_calc(in,od)
%
% This function is used to determine which individual moved to a new
% non-dominated solution. If so, the cognitive component used in the velocity
% equation is updated.
% Note: May want to try updated memory only when new solution dominates old
% solution in memory.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Jan. 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine New Pareto Solutions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if first iteration
if od.iter ~= 1

 % Initialize cognative influence parameters
  od.best(:,:,od.iter) = od.best(:,:,od.iter-1);
  od.best_fit(:,:,od.iter) = od.best_fit(:,:,od.iter-1);

 if 0
  % Mike Grant's:

  % Compare best fitness with current fitness
  diff_set = od.best_fit(:,:,od.iter) - od.fit(:,:,od.iter);

  % Determine characteristics of column in augmented matrix. Update cognitive
  % influcence parameters if new non-dominanated solution has been found. Note
  % this solution does not necessarily dominate the previous (may consider
  % changing).
  [r_big,c_big] = find(diff_set > 0);
  [r_small,c_small] = find(diff_set < 0);
  worse_only = setdiff(c_small,c_big);
  C = setdiff(1:1:in.num_ind,worse_only);

  % Update cognative influence parameters
  od.best(:,C,od.iter)     = od.pop(:,C,od.iter);
  od.best_fit(:,C,od.iter) = od.fit(:,C,od.iter);

 else

  % GFM
  % update best with the closest member in archive
  for i=1:in.num_ind
    dist = [norm(od.pop(:,i,od.iter)) 0];
    for j=find(~isnan(od.fit_arch(1,:,od.iter)))
      temp = norm(od.pop(:,i,od.iter)-od.pos_arch(:,j,od.iter));
      if temp < dist(1)
        dist = [temp j];
      end; 
    end; % FOR j

    % if we found a closest archive, use it 
    if dist(2) > 0
      od.best(:,i,od.iter) = od.pos_arch(:,dist(2),od.iter);
      od.best_fit(:,i,od.iter) = od.fit_arch(:,dist(2),od.iter);
    end;
  end; % FOR i
 end;
  
else
  
  % Initialize best fitness and best location parameters
  od.best(:,:,od.iter) = od.pop(:,:,od.iter);
  od.best_fit(:,:,od.iter) = od.fit(:,:,od.iter);
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [od] = determine_crowdedness(in,od)
%
% [od] = determine_crowdedness(in,od)
%
% This function is used to determine the crowdedness distance of the members of
% the archive.
%
% Input:
%   in - inputs from input deck [structure]
%   od - optimization variables of interest [structure]
%
% Output:
%   od - optimization variables of interest [structure]
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Feb. 2006


%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Determine size of archive
[R,C] = find(isnan(od.fit_arch(:,:,od.iter)) == 0);
C = unique(C);

% if archive is empty, return
if isempty(C)
  return;
end; 

% If two-objective optimization, then can determine endpoints of line. Endpoints
% are assigned an infinite crowdedness distance. This ensures the endpoints will
% always be picked when constructing the pool used by the social component of
% the velocity equation. This ensures the swarm will always try to expand the
% ends of the Pareto front.
if in.num_obj == 2

  % Determine endpoints of archive. Check to make sure consistent regardless of
  % dimension.
  for counter = 1 : 1 : in.num_obj

    [Y,I] = sort(od.fit_arch(counter,C,od.iter));
    pareto_lower(counter) = I(1);
    pareto_upper(counter) = I(end);

  end

  % Error check that Pareto curve always monotonically increasing in each
  % dimension.
  if sum(sort(pareto_lower)-sort(pareto_upper)) ~= 0
    
    error('Error: Pareto front not monotonically increasing! Algorithm invalid!');
    
  end
  
  % Determine if enough neighbors in archive
  if length(C) < 3
    
    % Account for not enough neighbors and if C is empty for average calcs below
    % (isempty is added to make sure don't divide by zero)
    num_ave = length(C) - 1;
    
  else
    
    % Use nearest 2 for crowdedness distance calculation
    num_ave = 2;
    
  end
  
elseif in.num_obj == 3
  
  % 3-D. Good luck in determining edges for Pareto front focus of expansion.
  % Determine if enough neighbors in archive
  if length(C) < 5
    
    % Account for not enough neighbors
    num_ave = length(C) - 1;
    
  else
    
    % Use nearest 4 for crowdedness distance calculation
    num_ave = 4;
    
  end
  
  % error('Not ready for 3-D yet!');
  
else
  
  error('Number of objectives can only be 2 or 3 in post_obj_calcs.m!');
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Crowdedness for Each Archive Member %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through each member of the archive
for counter = 1 : 1 : length(C)
  
  % Determine if boundary point of Pareto curve
  if (in.num_obj == 2) && (~isempty(find(pareto_lower == counter)))
    
    % Have boundary point on pareto front => assign infinite crowdedness
    % distance
    od.cd(1,C(counter),od.iter) = inf;
    
  else
    
    % Do not have boundary point 
    
    if num_ave == 0
      
      % Only one member of archive, no calculations necessary. Set crowdedness
      % distance to zero.
      od.cd(1,C(counter),od.iter) = 0;
      
    else

      % Other individuals in archive. Determine crowdedness distance.
    
      % Construct individual matrix to compare with archive
      individual = od.fit_arch(:,C(counter),od.iter);
      individual_mat = individual(:,ones(size(od.fit_arch(:,C,od.iter),2),1));

      % Compute relative positions
      rel_pos = od.fit_arch(:,C,od.iter) - individual_mat;

      % Determine spread in fitness space for each dimension
      spread = abs(max(rel_pos,[],2) - min(rel_pos,[],2));
      spread = spread(:,ones(size(rel_pos,2),1));

      % Normalize relative position values based on spread
      rel_pos = rel_pos./spread;

      % Sort normalized distances
      dist = sqrt(sum(rel_pos.^2,1));
      [dist_sorted,I] = sort(dist);

      % Determine crowdedness by averaging distance from two closest neighbors
      % in fitness space. Pareto front 2-D line.
      od.cd(1,C(counter),od.iter) = 1/num_ave*sum(dist_sorted(2:1+num_ave));
      
    end
      
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Average Crowdedness Distance %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute average crowdedness distance (not including Inf values)

% Determine set of valid crowdedness distance values
cd_set = od.cd(1,C,od.iter);

% Remove elements with values "Inf"
I = find(cd_set ~= inf);
cd_set = cd_set(I);

% Compute average crowdedness distance. Need 4 points to calculate mean. First 2
% are Inf since on boundaries, and need another 2 to calculate mean.
if length(cd_set) < 4

  % Not enough points to calculated average crowdedness distance. Set to 1.
  od.cd_mean(1,od.iter) = 1;
  
else

  % Enough points to calculate average crowdedness distance
  od.cd_mean(1,od.iter) = mean(cd_set);
  
end

return

