function out = particle_swarm
%
% particle_swarm utilizes the communication of all its particles in order to
%   determine a global maximum solution. The fitness of individual solutions is
%   specified in fitness_func.m. If the output is specified (i.e. there exists 
%   a calling function), then the output parameters are assigned:
%
% Output:
%
%   out -- Structure containing all of the swarm information
%     
%       out.ub -- Upper bound (see below)
%       out.lb -- Lower bound (see below)
%       out.solution_set -- History of the location of the particles in the 
%         swarm
%       out.best_fitness -- History of the best fitness in each iteration
%       out.num_iterations -- Number of iterations performed
%
% Particle Swarm Parameters:
%   
%   num_particles -- Determines the number of particles in the swarm. Increasing
%     the number of particles may help solution convergence but would slow
%     the processing of each iteration. A balance of swarm size and the
%     terminating criteria must be developed
%
%   num_neighbors -- Determines how many of the closest neighbors each
%     individual particle should communicate with
%
%   maxIterations -- Currently the terminating criteria of the swarm depends on
%     the number of iterations performed
%
%   deltaMin,deltaMax -- Constrains the minimum/maximum velocity of the 
%     particles in the search of the solution space
%
%   lb,ub -- Constrains the initial guess of the particles in the swarm
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine how function is used %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 0
  close all; clc;
  tic
end

%%%%%%%%%
% Input %
%%%%%%%%%

% Simlulation parameters
num_particles = 40;  % Number of particles in the swarm
num_neighbors = 10;  % Number of neighbors that will communicate with each particle
maxIterations = 100;  % Maximum number of iterations to perform
deltaMin = [-4.0 -4.0]';  % Lower bound of velocity
deltaMax = [4.0 4.0]';  % Upper bound of velocity
ub = [20 20]';  % Upper bound of trade-space
lb = [-20 -20]';  % Lower bound of trade-space

% Termination criteria, currently not in use (good for future work)
% Both variables determine what fitness 
delFitness = 0.1;  % Fitness-determined termination criteria
delIter = 10;  

% Constants
%% individuality and sociality
c1 = 2.0; % Maybe try 0.5 as well
c2 = 2.0; % Maybe try 0.5 as well
X = 1; % Constriction factor to limit velocity. Velocity is explicitly limited in code.
      % Hence, X is not needed. Setting X = 1 does not affect the velocity of the particles.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state',0); % Reset random number seed if want to have same initial conditions
iteration = 1;  % Iteration number
dim = length(deltaMin); % Determine number of dimensions of trade space

% Create a random initial guess grid that potentially covers the entire trade
% space
rand_matrix = rand(dim,num_particles);
diff_matrix = ub-lb;
current = rand_matrix.*diff_matrix(:,ones(size(rand_matrix,2),1))+lb(:,ones(size(rand_matrix,2),1));

% Create random initial velocities for all the particles in the swarm
rand_matrix = rand(dim,num_particles);
diff_matrix = deltaMax-deltaMin;
velocity = rand_matrix.*diff_matrix(:,ones(size(rand_matrix,2),1))+deltaMin(:,ones(size(rand_matrix,2),1));

next = current; % Next represents where the particles will move next. Initially next = current for simplicity.
neighbor = zeros(num_neighbors,num_particles); % Matrix used to keep track of neighbors.
best_fitness = 0*ones(1,maxIterations); % Vector used to store the best fitness of each generation.
bestSoFar = best_fitness(1)*ones(1,num_particles); % Vector used to store the best fitness each particle has seen.
best = zeros(dim,num_particles); % Matrix used to store the location of the best each particle has seen.

% Variables to save during the process of the swarm
best_save = zeros(dim,num_particles,maxIterations);
best_save(:,:,iteration) = best;
current_save = zeros(dim,num_particles,maxIterations);
current_save(:,:,iteration) = current;

% Set initial neighbors
for p = 1 : 1 : num_particles
  pcurrent = current(:,p);
  pcurrent = pcurrent(:,ones(size(current,2),1));
  dista = sqrt(sum((pcurrent-current).^2,1));
  [temp,I] = sort(dista);
  neighbor(:,p) = I(1:num_neighbors)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop of the Swarm %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while iteration <= maxIterations
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Move and Constrain Swarm %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Move particles to next location
  current = next;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Perform fitness operations for swarm %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Compute fitness of particles
  fitness = fitness_func(current);
  
  % Determine if particles moved to a higher fitness solution
  I = find(fitness-bestSoFar > 0);
  bestSoFar(I) = fitness(I);
  best(:,I) = current(:,I);
  best_fitness(iteration) = max(fitness);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Termination Criteria for Search %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  if abs(best_fitness(iteration)-best_fitness(iteration-1)) <= delFitness
%    X=particles(p).current(1)
%    Y=particles(p).current(2)
%    break;
%  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Determine new particle velocities %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Values constant per iteration
  r1 = rand;
  r2 = rand;
  w = 1.2 - iteration/maxIterations; % Inertia
  
  % Determine velocity of each particle
  for p = 1:num_particles
    
    % Set neighbors
    pcurrent = current(:,p);
    pcurrent = pcurrent(:,ones(size(current,2),1));
    dista = sqrt(sum((pcurrent-current).^2,1));
    [temp,I] = sort(dista);
    neighbor(:,p) = I(1:num_neighbors)';

    % Determine best neighbor of particle
    fit = fitness(neighbor(:,p));
    [temp,I] = sort(fit);
    n = neighbor(I(end),p);  % Index of best neighbor
    
    % Determine velocity of particle
    term1 = c1*r1*(best(:,p)-current(:,p));
    term2 = c2*r2*(best(:,n)-current(:,p));
    vel = X*(w*velocity(:,p) + term1 + term2);

    % Constrain velocity
    I = find(vel-deltaMin < 0);
    vel(I) = deltaMin(I);
    I = find(vel-deltaMax > 0);
    vel(I) = deltaMax(I);
      
    velocity(:,p) = vel;
    
    % Determine if particle would exceed bounds (approaches a wall).
    % If so, perform a perfectly elastic "bounce" and move particle exactly to
    % the wall.
    next(:,p) = current(:,p) + velocity(:,p);
    
    % Determine if upper bound exceeded. If so, place particle at wall.
    I = find(next(:,p) > ub);
    next(I,p) = ub(I);
    velocity(I,p) = -velocity(I,p); % Bounce
    
    % Determine if lower bound exceeded. If so, place particle at wall.
    I = find(next(:,p) < lb);
    next(I,p) = lb(I);
    velocity(I,p) = -velocity(I,p); % Bounce
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Save results and location of swarm %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  current_save(:,:,iteration) = current;
  best_save(:,:,iteration) = best;
  iteration = iteration + 1;
  
end

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

if nargout == 0
  % No calling function, just running genetic algorithm. Show output.
  
  % Data
  I = find(best_fitness==max(best_fitness));
  fprintf('\nBest Fitness Attained in Iteration %d out of %d',min(I),maxIterations);
  fprintf('\nBest Fitness Value: %g\n',best_fitness(min(I)));
  toc
  
  % Get fitness function data
  fitness_plot = fitness_func;

  % Plots
  figure(1); hold on;
  title('Position of Individuals in Population - Current');
  xlabel('x'); ylabel('y');
  %xlim([lb(1) ub(1)]);
  %ylim([lb(2) ub(2)]);
  [c,h] = contour(fitness_plot.X,fitness_plot.Y,fitness_plot.Z);

  for counter = 1 : 1 : length(current_save(1,1,:))
    h=plot(current_save(1,:,counter),current_save(2,:,counter),'.');
    pause(0.1);
    if counter ~= length(current_save(1,1,:))
      delete(h);
    end
  end
  
else
  % Calling function. Assign returning output.
  
  out.ub = ub;
  out.lb = lb;
  out.solution_set = current_save;
  out.best_fitness = best_fitness;
  out.num_iterations = maxIterations;
end
