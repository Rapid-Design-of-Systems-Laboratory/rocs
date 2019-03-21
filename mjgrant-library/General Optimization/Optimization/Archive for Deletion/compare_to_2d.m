function compare_to_2d

data = load('results.mat');
in = data.in;
od = data.od;

% Determine Pareto front (not just archive at end of optimization)
h = waitbar(0,'Generating Pareto Front...');
pareto_front = [];
for ctr = 1 : 1 : od.iter
  
  % Determine portion of current population that did not violate any constraints
  [R,C_fit] = find(isnan(od.fit(:,:,ctr)) == 0);
  C_fit = unique(C_fit); % Make sure no repetitions

  % Create augmented matrix combining fitness and Pareto front matrices
  augmented_matrix = [od.fit(:,C_fit,ctr) pareto_front];
  
  % Determine size of matrices
  [row_fit,col_fit] = size(od.fit(:,C_fit,ctr));
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
    fit = od.fit(:,C_fit(fitness_counter),ctr);
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
  
  pareto_keep = setdiff([1:1:size_pareto(2)],dominated_archive_set);

  % Remove dominated elements from Pareto front matrix
  pareto_front = pareto_front(:,pareto_keep);

  % Add Elements to Pareto Front
  pareto_front = [pareto_front od.fit(:,non_dominated_set,ctr)];
  
  waitbar(ctr/od.iter,h);

end
close(h);

plot3(pareto_front(1,:),pareto_front(2,:),pareto_front(3,:),'b.');

data = load('results_constrained_2d');
in = data.in;
od = data.od;

% Determine Pareto front (not just archive at end of optimization)
h = waitbar(0,'Generating Pareto Front...');
pareto_front = [];
for ctr = 1 : 1 : od.iter
  
  % Determine portion of current population that did not violate any constraints
  [R,C_fit] = find(isnan(od.fit(:,:,ctr)) == 0);
  C_fit = unique(C_fit); % Make sure no repetitions

  % Create augmented matrix combining fitness and Pareto front matrices
  augmented_matrix = [od.fit(:,C_fit,ctr) pareto_front];
  
  % Determine size of matrices
  [row_fit,col_fit] = size(od.fit(:,C_fit,ctr));
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
    fit = od.fit(:,C_fit(fitness_counter),ctr);
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
  
  pareto_keep = setdiff([1:1:size_pareto(2)],dominated_archive_set);

  % Remove dominated elements from Pareto front matrix
  pareto_front = pareto_front(:,pareto_keep);

  % Add Elements to Pareto Front
  pareto_front = [pareto_front od.fit(:,non_dominated_set,ctr)];
  
  waitbar(ctr/od.iter,h);

end
close(h);

view(0,90);
hold on;
plot(pareto_front(1,:),pareto_front(2,:),'r.');
legend('3D Pareto Surface Projection','2D Pareto Front');

return

