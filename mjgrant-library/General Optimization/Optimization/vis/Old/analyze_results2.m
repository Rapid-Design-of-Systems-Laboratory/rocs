function analyze_results

%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Determine if create movie file
make_movie = false;

% Determine if pausing in real-time figure updating
pausing = true;

% Initialize plot information
pos1 = [0.05 0.50 0.46 0.45];
pos2 = [0.08 0.06 0.33 0.33];
pos3 = [0.58 0.50 0.40 0.45];
pos4 = [0.50 0.08 0.33 0.33];
pos5 = [0.83 0.08 0.15 0.30];

% Set figure
figure;
set(gcf,'color','white');
set(gcf,'Position',[1 31 1023 663]);

%%%%%%%%%%%%%%%%%
%% Import Data %%
%%%%%%%%%%%%%%%%%

load('results.mat');

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

% Determine convergence
if od.iter == in.max_iter

    fprintf('\nSolution not necessarily converged.\n');
    fprintf('Maximum number of iterations performed (%g iterations)\n',od.iter);
  
  else
  
    fprintf('\nSolution Converged\n');
    fprintf('Number of Iterations: %g\n',od.iter);
  
end

% Execution time
fprintf('Time elapsed: %g\n',od.time);

% Determine number of objectives in optimization
if in.num_obj == 1

  % Single-objective optimization

  % Negate fitness if performing maximization
  if in.maximization
    od.best_fit = -od.best_fit;
    od.fit = -od.fit;
  end

  [Y,I] = min(od.best_fit);
  I = min(I);
  [Y,p] = min(od.fit(1,:,I));
  p = min(p);

  if in.maximization
    od.best_fit = -od.best_fit;
    od.fit = -od.fit;
  end

  % Data
  fprintf('\nOptimal Solution Attained in Iteration %d out of %d',I,od.iter);
  fprintf('\nOptimal Solution Value: %g',od.best_fit(I));
  fprintf('\nTrade Space Location:\n');
  fprintf('      %g\n',od.pop(:,p,I));
  
else
  
  % Multi-objective optimization
  
  % Check to make sure more than one non-dominated solution
  if length(find(isnan(od.cd) == 0)) == 1
    
    % Only one non-dominated solution. Objectives may not be in competition
    fprintf('Warning! Only one non-dominated solution found!\n');
    
  end

end

fprintf('\n\n');

%%%%%%%%%%%
%% Plots %%
%%%%%%%%%%%

% Determine number of objectives in optimization
if in.num_obj == 1
  
  % Single-objective optimization - generic plots
  
%   % Contour of positions
%   x = [];
%   y = [];
%   fitness = [];
% 
%   if od.dim ~= 1
%   
%     % Construct fitness contour
%     for iter_counter = 1 : 1 : od.iter
%       fitness = [fitness od.fit(:,:,iter_counter)];
%       x = [x od.pop(1,:,iter_counter)];
%       y = [y od.pop(2,:,iter_counter)];
%     end
% 
%     x_min = min(x);
%     x_max = max(x);
%     y_min = min(y);
%     y_max = max(y);
% 
%     step_x = (x_max - x_min)/100;
%     step_y = (y_max - y_min)/100;
%     xg = [x_min : step_x : x_max];
%     yg = [y_min : step_y : y_max];
% 
%     [X,Y] = meshgrid(xg,yg);
%     warning off; % Turn off warning for duplicate data
%     fitness_plot = griddata(x,y,fitness,X,Y);
%     warning on;
%     
%   end
% 
%   % Create subplot of current solution state
%   subplot('position',pos1); hold on; set(gca,'box','on'); grid on;
%   title('Position of Individuals in Population - Current');
%   xlabel('x');
%   xlim([in.pos_min(1) in.pos_max(1)]);
%   if od.dim ~= 1
%     ylabel('y');
%     ylim([in.pos_min(2) in.pos_max(2)]);
%     contour(X,Y,fitness_plot);
%   end
  
  if strcmp(in.method,'pso')
    % PSO - Create subplot of inertia
    subplot('position',pos2); set(gca,'box','on'); grid on; hold on;
    title('Inertia vs. Iteration');
    xlabel('Iteration'); ylabel('Inertia');
    plot([1:1:od.iter],od.w(1:od.iter));
    axis_info = get(gca);
  elseif strcmp(in.method,'ga')
    % GA - Do nothing
  end

%   % Create subplot of solution state history
%   subplot('position',pos3); set(gca,'box','on'); hold on; grid on;
%   title('Position of Individuals in Population - History');
%   xlabel('x');
%   xlim([in.pos_min(1) in.pos_max(1)]);
%   if od.dim ~= 1
%     ylabel('y');
%     ylim([in.pos_min(2) in.pos_max(2)]);
%     contour(X,Y,fitness_plot);
%   end

  % Create subplot of best fitness curve
  subplot('position',pos4); set(gca,'box','on'); grid on; hold on;
  title('Best Fitness vs. Iteration');
  xlabel('Iteration'); ylabel('Best Fitness');
  plot([1:1:od.iter],od.best_fit(1:od.iter));
  axis_info2 = get(gca);

  % Create "subplot" of text information
  subplot('position',pos5); set(gca,'Visible','off'); hold on;

  % Create movie object and start recording
  if make_movie
    aviobj = avifile('optimization_movie.avi','Compression','Indeo5','FPS',8);
%       aviobj = avifile(movie(optim_counter).name,'FPS',8);
  end
  
  % Plot data and create movie frame by frame
  for counter = 1 : 1 : od.iter
    % Plot particle swarm info
    subplot('position',pos1);
    if od.dim == 1
      h1 = plot(od.pop(1,:,counter),zeros(1,length(od.pop(1,:,counter))),'.');
    else
      h1 = plot(od.pop(1,:,counter),od.pop(2,:,counter),'.');
    end
    subplot('position',pos4);
    h4 = plot([counter counter],axis_info2.YLim,'r');
    subplot('position',pos3);
    if od.dim == 1
      h2 = plot(od.pop(1,:,counter),zeros(1,length(od.pop(1,:,counter))),'.');
    else
      h2 = plot(od.pop(1,:,counter),od.pop(2,:,counter),'.');
    end
    subplot('position',pos5);
    h5 = text(0.1,0.8,['Iteration: ',num2str(counter)]);
    if strcmp(in.method,'pso')
      subplot('position',pos2);
      h6 = plot([counter counter],axis_info.YLim,'r');
    elseif strcmp(in.method,'ga')
      h6 = [];
    end
    if make_movie
      h = gcf;
      aviobj = addframe(aviobj,h); % Add figure to movie
    end
    if pausing; pause(0.1); end
    if counter ~= od.iter
      delete(h1,h4,h5,h6);
    end
  end
  
  if make_movie
    aviobj = close(aviobj);
  end

  if strcmp(in.method,'pso')
  
    % Fitness Coefficient of Variance vs. Iteration
    figure;
    plot([1:od.iter],od.cov(1:od.iter));
    title('Fitness COV vs. Iteration');
    xlabel('Iteration');
    ylabel('Coefficient of Variance');
    grid on;

    % Velocity vs. Iteration - Each Dimension
    for dim_counter = 1 : 1 : od.dim
      figure; hold on; grid on;
      % Formatting
      title(['Velocity of Each Individual - Dim ',num2str(dim_counter)]);
      xlabel('Iteration');
      ylabel('Velocity');
      % Plot velocity of each individual for each iteration
      for ind_counter = 1 : 1 : in.num_ind
        data(1:od.iter) = od.vel(dim_counter,ind_counter,[1:od.iter]);
        plot([1:od.iter],data,'.');
      end
%       % Size window to bounds
%       axis([1 od.iter in.vel_min(dim_counter) in.vel_max(dim_counter)]);
    end
    
  end
  
else
  
  % Multi-objective particle swarm

  % Trade space plot
  subplot('position',pos1); set(gca,'box','on');
  
  % Inertia plot
  subplot('position',pos2); hold on; set(gca,'box','on');
  title('Inertia vs. Iteration');
  xlabel('Iteration'); ylabel('Inertia');
  plot([1:1:od.iter],od.w(1:od.iter));
  axis_info = get(gca);

  % Fitness space plot
  subplot('position',pos3); set(gca,'box','on');

  % Create subplot of average crowdedness distance
  subplot('position',pos4); set(gca,'box','on'); grid on; hold on;
  title('Mean Normalized Crowdedness Distance vs. Iteration');
  xlabel('Iteration'); ylabel('Mean Normalized Crowdedness Distance');
  plot([1:1:od.iter],od.cd_mean(1:od.iter));
  axis_info2 = get(gca);

  % Text box
  subplot('position',pos5); set(gca,'Visible','off'); hold on;

  % Create movie object and start recording
  if make_movie
    aviobj = avifile('optimization_movie.avi','Compression','Indeo5','FPS',8);
%       aviobj = avifile(movie(optim_counter).name,'FPS',8);
  end
  
  for counter = 1 : 1 : od.iter

    % Plot front in trade space
    subplot('position',pos1);
    title('Evolution of Pareto Front and Swarm - Trade Space');
    xlabel('X1'); ylabel('X2');
    if od.dim == 2
      % Set plot size
      axis([in.pos_min(1) in.pos_max(1) in.pos_min(2) in.pos_max(2)]);
      % Plot front
      h1 = plot(od.pos_arch(1,:,counter),od.pos_arch(2,:,counter),'b.');
      hold on; grid on;
      % Plot swarm
      h2 = plot(od.pop(1,:,counter),od.pop(2,:,counter),'rx');
    elseif od.dim == 3
      zlabel('X3');
      % Set plot size
      axis([in.pos_min(1) in.pos_max(1) in.pos_min(2) in.pos_max(2) ...
      in.pos_min(3) in.pos_max(3)]);
      % Plot front
      h1 = plot3(od.pos_arch(1,:,counter),od.pos_arch(2,:,counter),...
        od.pos_arch(3,:,counter),'b.');
      hold on; grid on;
      % Plot swarm
      h2 = plot3(od.pop(1,:,counter),od.pop(2,:,counter),...
        od.pop(3,:,counter),'rx');
    end

    % Plot front in fitness space
    subplot('position',pos3);
    title('Evolution of Pareto Front and Swarm - Fitness Space');
    xlabel('J_1'); ylabel('J_2');
    if in.num_obj == 2
      % Set plot size
      trade_min = min(od.fit_arch(:,:,counter),[],2);
      trade_max = max(od.fit_arch(:,:,counter),[],2);
      axis([trade_min(1) trade_max(1) trade_min(2) trade_max(2)]);
      % Plot front
      h3 = plot(od.fit_arch(1,:,counter),od.fit_arch(2,:,counter),'b.');
      hold on; grid on;
      % Plot swarm
      h4 = plot(od.fit(1,:,counter),od.fit(2,:,counter),'rx');
    elseif in.num_obj == 3
      zlabel('J_3');
      % Set plot size
      trade_min = min(od.fit_arch(:,:,counter),[],2);
      trade_max = max(od.fit_arch(:,:,counter),[],2);
      axis([trade_min(1) trade_max(1) trade_min(2) trade_max(2) ...
        trade_min(3) trade_max(3)]);
      % Plot front
      h3 = plot3(od.fit_arch(1,:,counter),od.fit_arch(2,:,counter),...
        od.fit_arch(3,:,counter),'b.');
      hold on; grid on;
      % Plot swarm
      h4 = plot3(od.fit(1,:,counter),od.fit(2,:,counter),...
        od.fit(3,:,counter),'rx');
    end

    % Plot markers in inertia and crowdedness distance plots
    subplot('position',pos2);
    h5 = plot([counter counter],axis_info.YLim,'r');
    subplot('position',pos4);
    h6 = plot([counter counter],axis_info2.YLim,'r');

    % Text information
    subplot('position',pos5);
    h7 = text(0.1,0.8,['Iteration: ',num2str(counter)]);

    if make_movie
      h = gcf;
      aviobj = addframe(aviobj,h); % Add figure to movie
    end
    if pausing; pause(0.1); end
    % Pause and delete, if necessary
    if counter ~= od.iter
      delete(h1,h2,h3,h4,h5,h6,h7);
    end

  end
  
  if make_movie
    aviobj = close(aviobj);
  end
  
  % Plot parallel coordinate system
  
  [m,n] = size(od.fit_arch(:,:,od.iter));
  max_vals = max(od.fit_arch(:,:,od.iter),[],2);
  max_vals = ones(m,1);
  arch_norm = od.fit_arch(:,:,od.iter) ./ ...
    max_vals(:,ones(size(od.fit_arch(:,:,od.iter),2),1));

  figure; hold on;
  set(gcf,'color','white');
  set(gcf,'Position',[1 31 1023 663]);
  
  for member = 1 : 1 : n
    
    plot([1:1:m],arch_norm(:,member));

  end
  set(gca,'xtick',[0:1:in.num_obj]) 
  xlabel('Objective');
  ylabel('J');
  title('Parallel Coordinate System of Pareto Front');
  
  % Big plot of Pareto Front
  figure;
  set(gcf,'color','white');
  set(gcf,'Position',[1 31 1023 663]);
  if in.num_obj == 2
    plot(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter),'b.');
  elseif in.num_obj == 3
    plot3(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter), ...
      od.fit_arch(3,:,od.iter),'b.');
    zlabel('J3');
  end
  title('Pareto Front');
  xlabel('J1');
  ylabel('J2');
  grid on;
  
  
  
end

