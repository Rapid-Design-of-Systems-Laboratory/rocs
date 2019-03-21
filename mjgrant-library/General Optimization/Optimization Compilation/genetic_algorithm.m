function out = genetic_algorithm
%
% genetic_algorithm utilizes the evolution of all its members in order to
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
%       out.solution_set -- History of the location of the population
%       out.best_fitness -- History of the best fitness in each generation
%       out.num_iterations -- Number of iterations performed
%
% Genetic Algorithm Parameters:
%   
%   pop_size -- Determines the number of individuals in the generation. Increasing
%     the number of individuals may help solution convergence but would slow
%     the processing of each generation. A balance of generation size and the
%     terminating criteria must be developed
%
%   max_generation -- Currently the terminating criteria of the GA depends on
%     the number of iterations performed
%
%   lowtrait,hightrait -- Constrains the minimum/maximum admissible values of
%     the independent variables of the GA in the search of the solution space
%
%   elitism -- Determines if the best member of a generation should be
%     transferred directly into the next generation
%  
%   cross_prob -- Determines the probability that a parent will mate with
%     another individual in the population to form children for the next
%     generation. If the parent does not mate, then the parent is directly
%     inserted into the next generation.
%
%   mutate_prob -- Determines the probability that any given trait will mutate.
%     The mutation, by design of the code, will be bounded in the admissible trade
%     space.
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

pop_size = 40;  % Population size - constant for each generation
hightrait = [20 20]';  % Upper bound for traits
lowtrait = [-20 -20]'; % Lower bound for traits
max_generation = 100; % Maximum number of generations to run
elitism = 1;  % Elitism on/off = 1/0
cross_prob = 0.8;  % Probability of crossover
mutate_prob = 0.05;  % Probability of mutation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_traits = length(hightrait); % Pick a matrix for number of traits
rand('state',0); % Reset random number seed if want to have same initial conditions
pop = zeros(num_traits,pop_size,max_generation); % Population matrix
fitness = zeros(1,pop_size);  % Fitness values for one generation
bestfitness = zeros(1,max_generation); % Best fitness values for all generations
popcount = 1; % Initialize population count

% Construct random initial guess for GA. The random guess for the initial 
% population should have the potential of covering the entire trade space 
% within the bounds specified by the input.
difftrait = (hightrait-lowtrait);
randmat = rand(num_traits,pop_size,popcount);
current_pop = randmat.*difftrait(:,ones(size(randmat,2),1))+ ...
    lowtrait(:,ones(size(randmat,2),1));

%%%%%%%%%%%%%%%%%%%%%%%
% Main loop of the GA %
%%%%%%%%%%%%%%%%%%%%%%%

while popcount <= max_generation

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Constrain traits to valid search space %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % High saturation limit exceeded
  [R,C] = find(current_pop > hightrait(:,ones(size(current_pop,2),1)));
  current_pop(R,C) = hightrait(R);
  
  % Low saturation limit exceeded
  [R,C] = find(current_pop < lowtrait(:,ones(size(current_pop,2),1)));
  current_pop(R,C) = lowtrait(R);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Perform fitness operations for generation %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Determine fitness of each individual in generation
  fitness = fitness_func(current_pop);
  
  % Determine the most fit individual
  [bestfitness(popcount),bestmember]=max(fitness);
  
  % Make all fitness values positive and start at zero for roulette selection
  fitness_abs = fitness - min(fitness);
  sumfitness = sum(fitness_abs);

  %%%%%%%%%%%%%%%%%%%%%
  % Determine parents %
  %%%%%%%%%%%%%%%%%%%%%
  
  % Form the mating pool by randomly selecting individuals to mate using a
  % roulette wheel approach. Individuals with high performance function values 
  % are more likely to mate. If elitism is on, the best individual is directly 
  % copied to the next generation and must be assigned as a parent.
  
  for pop_member = 1:pop_size
      
	if elitism ==1 && pop_member==bestmember  
      % If elitism on and have the elite member, copy elite member
      % into the next generation.
      parent_chrom(:,pop_member)=current_pop(:,pop_member);
	else
      pointer=rand*sumfitness; % This spins the roulette wheel      
      member_count=1;        % Initialization
      total=fitness_abs(1);
      % Determine which member pointer is identifying
      while total < pointer,                
         member_count=member_count+1;
         total=total+fitness_abs(member_count);
      end
      % Make the parent chromosome
      parent_chrom(:,pop_member)=current_pop(:,member_count);
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Determine children - Perform crossover %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Reproduction is performed using crossover. The point for crossover is
  % randomly chosen for each mating pair.
  
  for parent_number1 = 1:pop_size,  % parent_number1 is individual who gets to mate
      
    % Determine if elitism is on
	if elitism ==1 && parent_number1==bestmember 
      % Have the elite member - direct copy
	  child(:,parent_number1)=parent_chrom(:,parent_number1);
    else
      % Randomly determine mate from parent set
	  parent_number2=parent_number1;  % Initialize who the mate is
      % Iterate until find a random mate other than itself
	  while parent_number2 == parent_number1   
        parent_number2 = floor(rand*pop_size)+1; 
      end
      % Determine if crossover should be performed
      if cross_prob > rand           
        % Perform crossover by randomly choosing  site for crossover and 
        % make it a valid integer number for a site
        site = floor(rand*num_traits)+1;
        % Form the child by the swapping of genetic material between the parents
        child(1:site,parent_number1)=parent_chrom(1:site,parent_number1);
        child(site+1:num_traits,parent_number1)=...
          parent_chrom(site+1:num_traits,parent_number2);
      else
        % No crossover occurs
        % Copy non-crossovered chromosomes into next generation
        % In this case, simply take one parent and make them the child.
        child(:,parent_number1)=parent_chrom(:,parent_number1);
      end
    end  % End the "if ELITISM..." statement
    
  end  % End "for parent_number1=..." loop
  
  %%%%%%%%%%%%%%%%%%%
  % Mutate children %
  %%%%%%%%%%%%%%%%%%%
  
  % Randomly mutate children.  Used to enhance global searching capability.
  
  for pop_member= 1:pop_size,

    % Determine if elitism is on
    if elitism ==1 && pop_member==bestmember  
      % Have the elite member, DO NOT mutate the elite member
	  child(:,pop_member)=child(:,pop_member);
    else
      % Randomly mutate traits of children
	  for site = 1:num_traits

        if mutate_prob > rand
          % Mutate trait
          rand_gene=rand*(hightrait(site)-lowtrait(site))+lowtrait(site); % Create a random trait
          % If it is the same as the one already there then generate another random trait
          while child(site,pop_member) == rand_gene
            rand_gene=rand*(hightrait(site)-lowtrait(site))+lowtrait(site); % Create a random trait
          end;
          % If trait is not the same, then mutate
          child(site,pop_member)=rand_gene;
            
        end  % End "if MUTAT_PROB > rand ... 
	  end  % End for site... loop
	end  % End "if ELITISM..."
  end  % End for pop_member loop
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Create the next generation %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  current_pop=child;
  pop(:,:,popcount)=child;      % Create next generation (children become parents)
  popcount=popcount+1;    % Increment to the next generation
  
end  % End "for popcount=..." loop - the main loop.

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

if nargout == 0
  % No calling function, just running genetic algorithm. Show output.

  % Data
  I = find(bestfitness==max(bestfitness));
  fprintf('\nBest Fitness Attained in Generation %d out of %d',min(I),max_generation);
  fprintf('\nBest Fitness Value: %g\n',bestfitness(min(I)));
  toc

  % Get fitness function data
  fitness_plot = fitness_func;

  % Plots
  figure(1); hold on;
  title('Position of Individuals in Population - Current');
  xlabel('x'); ylabel('y');
  xlim([lowtrait(1) hightrait(1)]);
  ylim([lowtrait(2) hightrait(2)]);
  [c,h] = contour(fitness_plot.X,fitness_plot.Y,fitness_plot.Z);
  for counter = 1 : 1 : length(pop(1,1,:))
    h=plot(pop(1,:,counter),pop(2,:,counter),'.');
    pause(0.1);
    if counter ~= length(pop(1,1,:))
       delete(h);
    end
  end

else
  % Calling function. Assign returning output.
  
  out.solution_set = pop;
  out.best_fitness = bestfitness;
  out.lb = lowtrait;
  out.ub = hightrait;
  out.num_iterations = max_generation;
end




