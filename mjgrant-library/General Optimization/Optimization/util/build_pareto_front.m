function pareto_front = build_pareto_front(fit_set,wait_str,maximization)
%
% pareto_front = build_pareto_front(fit_set)
%
% This function is used to build the Pareto front from the objective information
% provided in fit_set. The Pareto front information is stored in pareto_front.
%

% Author NASA JSC/DM42 - Michael J. Grant in July 2006


% Determine how many iterations were performed in optimization
size_info = size(fit_set);
num_iter = size_info(end);

% Initialize
if ~exist('wait_str','var') || isempty(wait_str)
  wait_str = 'Building the Pareto Front...';
end

% Turn maximization objectives into minimization
if ~exist('maximization','var') | isempty(maximization)
  maximization = false*ones(1,size_info(1));
end
I_max = find(maximization == true);
fit_set(I_max,:,:) = -fit_set(I_max,:,:);
  
h = waitbar(0,wait_str);
pareto_front = [];

% Loop through each iteration to build the Pareto front
for ctr = 1 : 1 : num_iter
  
  % Determine portion of current population that did not violate any constraints
  [R,C_fit] = find(isnan(fit_set(:,:,ctr)) == 0);
  C_fit = unique(C_fit); % Make sure no repetitions

  % Create augmented matrix combining fitness and Pareto front matrices
  augmented_matrix = [fit_set(:,C_fit,ctr) pareto_front];
  
  % Determine size of matrices
  [row_fit,col_fit] = size(fit_set(:,C_fit,ctr));
  [row_augmented,col_augmented] = size(augmented_matrix);

  % Vector of individuals that are non-dominated
  non_dominated_set = [];

  % Vector of members in augmented matrix that are dominated by individual
  dominates = []; 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Test for Pareto Optimality %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Loop through each individual to test dominance
  for fitness_counter = 1 : 1 : col_fit

    % Create matrix repeating fitness values for individual (for matrix math)
    fit = fit_set(:,C_fit(fitness_counter),ctr);
    fit_new = fit(:,ones(size(augmented_matrix,2),1));

    % Create matrix of relative fitness values
    diff_set = augmented_matrix - fit_new;

    % Determine characteristics of column in augmented matrix (for dominance)
    [r_big,c_big] = find(diff_set > 0);
    [r_small,c_small] = find(diff_set < 0);
    c_small = unique(c_small);
    c_big = unique(c_big);

    % Determine dominance
    if length(c_big) == (col_augmented - 1)
    
      % Individual is not dominated by another other individual. Determine which 
      % members of the augmented matrix are dominated. Dominates includes values
      % from entire augmented matrix.
      dominates = [dominates setdiff(c_big,c_small)'];
    
      % Store individual as non-dominated solution
      non_dominated_set = [non_dominated_set C_fit(fitness_counter)];
    
    end

  end

  % Dominates includes values from entire augmented matrix. Trim to determine
  % those solutions from Pareto front matrix that are dominated.
  dominates = dominates - col_fit;
  I = find(dominates > 0);
  size_pareto = size(pareto_front);
  C = 1 : 1 : size_pareto(2);
  dominated_archive_set = C(unique(dominates(I)));
  
  % Determine which elements of the Pareto front to keep
  pareto_keep = setdiff([1:1:size_pareto(2)],dominated_archive_set);

  % Remove dominated elements from Pareto front matrix
  pareto_front = pareto_front(:,pareto_keep);

  % Add Elements to Pareto Front
  pareto_front = [pareto_front fit_set(:,non_dominated_set,ctr)];
  
  % Update waitbar
  waitbar(ctr/num_iter,h);

end

% Restore values if maximization performed
fit_set(I_max,:,:) = -fit_set(I_max,:,:);
pareto_front(I_max,:,:) = -pareto_front(I_max,:,:);

% Close waitbar
close(h);

return

