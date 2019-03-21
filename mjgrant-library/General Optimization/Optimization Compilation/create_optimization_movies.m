function create_optimization_movies

clear all; close all; clc;

% Inputs -----------------------------------------------------------------------
% Determine what optimization techniques are to be used
%   Genetic Algorithm (GA) = 1
%   Particle Swarm (Swarm) = 2
optimization_set = [1];
make_movie = true;
pausing = false;

% Set GA info
graphing(1).current = 'Position of Individuals in Population - Current';
graphing(1).history = 'Position of Individuals in Population - History';
movie(1).name = 'genetic_algorithm.avi';
% Set particle swarm info
graphing(2).current = 'Position of Particles in Swarm - Current';
graphing(2).history = 'Position of Particles in Swarm - History';
movie(2).name = 'particle_swarm.avi';

% Initialize plot information
pos1 = [0.05 0.50 0.46 0.45];
pos2 = [0.08 0.06 0.33 0.33];
pos3 = [0.58 0.50 0.40 0.45];
pos4 = [0.50 0.08 0.33 0.33];
pos5 = [0.83 0.08 0.15 0.30];

% Create Movie with Plots ------------------------------------------------------

% Get Fitness Function Data
fitness_plot = fitness_func;

for optim_counter = optimization_set
    
    % Determine which optimization technique to run
    switch optim_counter
      case 1
        optim_info = genetic_algorithm; % Run GA
      case 2
        optim_info = particle_swarm; % Run swarm
    end
    
    % Set figure
    figure(optim_counter);
    set(gcf,'color','white');
    set(gcf,'Position',[1 31 1023 663]);
    
    % Create subplot of current solution state
    subplot('position',pos1); hold on; set(gca,'box','on')
    title(graphing(optim_counter).current);
    xlabel('x'); ylabel('y');
    xlim([optim_info.lb(1) optim_info.ub(1)]);
    ylim([optim_info.lb(2) optim_info.ub(2)]);
    x = [optim_info.lb(1) : .5 : optim_info.ub(1)];
    y = [optim_info.lb(2) : .5 : optim_info.ub(2)];
    [c,h] = contour(fitness_plot.X,fitness_plot.Y,fitness_plot.Z); colorbar;
    data = get(gca);
    
    % Create subplot of 3-D fitness function
    subplot('position',pos2); set(gca,'box','on'); grid on;
    surf(fitness_plot.X,fitness_plot.Y,fitness_plot.Z);
    title('Fitness Function');
    xlabel('x'); ylabel('y'); zlabel('Fitness');
    
    % Create subplot of solution state history
    subplot('position',pos3); set(gca,'box','on'); hold on;
    title(graphing(optim_counter).history);
    xlabel('x'); ylabel('y');
    xlim([optim_info.lb(1) optim_info.ub(1)]);
    ylim([optim_info.lb(2) optim_info.ub(2)]);
    [c,h] = contour(fitness_plot.X,fitness_plot.Y,fitness_plot.Z);
    
    % Create subplot of fitness curve
    subplot('position',pos4); set(gca,'box','on'); grid on; hold on;
    title('Best Fitness vs. Iteration');
    xlabel('Iteration'); ylabel('Fitness');
    plot([1:1:optim_info.num_iterations],optim_info.best_fitness);
    axis_info2 = get(gca);
    
    % Create "subplot" of text information
    subplot('position',pos5); set(gca,'Visible','off'); hold on;
    
    % Create movie object and start recording
    if make_movie
      aviobj = avifile(movie(optim_counter).name,'Compression','Indeo5','FPS',8);
%       aviobj = avifile(movie(optim_counter).name,'FPS',8);
    end
    
    % Plot data and create movie frame by frame
    for counter = 1 : 1 : optim_info.num_iterations
      % Plot genetic algorithm info
      subplot('position',data.Position);
      h1=plot(optim_info.solution_set(1,:,counter),optim_info.solution_set(2,:,counter),'.');
      subplot('position',pos4);
      h4 = plot([counter counter],axis_info2.YLim,'r');
      subplot('position',pos3);
      h2 = plot(optim_info.solution_set(1,:,counter),optim_info.solution_set(2,:,counter),'.');
      subplot('position',pos5);
      h5 = text(0.1,0.8,['Iteration: ',num2str(counter)]);
      if make_movie
        h = gcf;
        aviobj = addframe(aviobj,h); % Add figure to movie
      end
      if pausing; pause(0.1); end
      if counter ~= optim_info.num_iterations
        delete(h1,h4,h5);
      end
    end
    
    if make_movie
      aviobj = close(aviobj);
    end

end
