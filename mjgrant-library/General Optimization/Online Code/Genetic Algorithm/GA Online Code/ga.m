%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Genetic Algorithm (GA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kevin Passino
%   Version: 7/21/98, latest change 10/22/03 with correction from 
%    Phan Tran Ho Truc, Mr., B. Eng., lecturer in HCM City University of Technology
%    to lines 379-384 of old program (lines 380-381 of this program)
%
% Notes: This program has evolved (hopefully to achieve a higher fitness!)
% over time using programming ideas from several persons including 
% LaMoyne Porter, Will Lennon, Jonathan Cook, and Jim Gremling.  
% 
% This program simulates the optimzation of a simple function with
% a base-10 genetic algorithm.  First, we outline some background 
% information on GAs and define some terminology.  Following this,
% we explain how to use the program for your own purposes and as an
% example we show how to solve an optimization problem.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Background Information on GAs:
%
%    A genetic algorithm is an optimization method that
% manipulates a string of numbers in a manner similar to how chromosomes  
% are changed in biological evolution.  An initial population made up of 
% strings of numbers is chosen at random or is specified by the 
% user. Each string of numbers is called a "chromosome" or an 
% "individual," and each number slot is called a "gene." A 
% set of chromosomes forms a population.  Each chromosome 
% represents a given number of traits which are the actual 
% parameters that are being varied to optimize the "fitness function".  
% The fitness function is a performance index that we seek to maximize.  
%
%    The operation of the GA proceeds in steps. Beginning with the initial 
% population, "selection" is used to choose which chromosomes 
% should survive to form a "mating pool."  Chromosomes are chosen based on 
% how fit they are (as computed by the fitness function) relative
% to the other members of the population.  More fit individuals end up 
% with more copies of themselves in the mating pool so that they
% will more significantly effect the formation of the next 
% generation.  Next, several operations are taken on the mating
% pool.  First, "crossover" (which represents mating, the exchange
% of genetic material) occurs between parents.
% To perform crossover, a random spot is picked 
% in the chromosome, and the genes after this spot are switched 
% with the corresponding genes of the other parent.  Following this, 
% "mutation" occurs.  This is where some genes are randomly changed 
% to other values.  After the crossover and mutation operations
% occur, the resulting strings form the next generation
% and the process is repeated.  A termination criterion
% is used to specify when the GA should end (e.g., the maximum 
% number of generations or until the fitness stops increasing).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Some GA Definitions for this Program:
%  
% A GENE is a digit position that can take on certain values
%   Ex: In a base-10 algorithm a gene can hold any number between 0 and 9.
%
% A CHROMOSOME is a string of genes.
%   Ex: In base-10 1234567890 could be a chromosome of length 10.
%
% A TRAIT is a decimal number which is decoded from a chromosome.
%   Normally, a chromosome is a concatenation of several TRAITS.
% 
% An INDIVIDUAL is the object that the GA is attempting to optimize.
%  An individual is described by its chromosome.
%  The individual's traits determine its fitness.
%
% A POPULATION is a set of individuals (set of chromosomes).
% 
% FITNESS FUNCTION is the objective function to be optimized which provides
%  the mechanism for evaluating each string (we maximize the fitness 
%  function).
% 
% SELECTION: a fitter string receives a higher number of offspring 
%  and thus has a higher chance of surviving in the subsequent generation.
%
% CROSSOVER is an operation which swaps genes between two chromosome.
% 
% MUTATION is flipping a digit in the chromosome after crossover.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   How to use the program for your own purposes:
%
% The main thing that needs to be changed in ga.m is the fitness function 
% and a few parameters.  It is now set to minimize a function z=f(x,y)
% that is a sum of scaled translated Gaussian distributions with x and y 
% between 0 and 30 (specified by the elements of HIGHTRAIT and LOWTRAIT).
% The function to be optimized must be changed in the section specified 
% in the program below.
%
%    In addition, the following parameters should be modified according
% to your needs:
%
% 1. NUM_TRAITS - the number of traits represented by each chromosome
%                  (e.g. if we were trying to optimize the above function
%				   we could have NUM_TRAITS=2).
%                   Example: Chromosome: 123456789012
%                        Possible two traits: 123456, 789012
%
% 2. HIGHTRAIT,LOWTRAIT -  arrays with NUM_TRAITS elements with each element
%                 specifying the upper (lower) limit on each trait (often
% 			      known a priori for an optimization problem).
%                 For example, we may only want to find the optimum values
%				  of a function over a certain region.  For instance, 
%				  we could choose HIGHTRAIT=[30 30] and LOWTRAIT=[0 0]
%				  for the above optimization problem.
%
% 3. SIG_FIGS - an array with NUM_TRAITS elements with each element 
%                specifying the number of significant figures in each 
%                trait (there must be at least one significant figure).
%                Example: SIG_FIGS=[6 6] with above possible traits.
%
% 4. DECIMAL - an array with NUM_TRAITS elements with each element
%               specifying the number of digits to the left of the decimal
%               point for each trait.  Choose DECIMAL to be the 
%               largest number integer part in HIGHTRAIT or LOWTRAIT
%               to avoid saturation.
%               (Example: HIGHTRAIT=[12.9 9.9], LOWTRAIT=[5.5 -5.5],
%                    DECIMAL=[2 1] so for our optimization problem 
%                    let DECIMAL=[2 2])
%
%   With this we can explain how to "decode" a chromosome:
%
%    The first gene of each trait determines the sign of the
% trait.  For our base-10 algorithm if this gene is 0-4, the trait is
% negative and if this gene is 5-9, the trait is positive.  
% The remaining genes determine the size of the trait.  The 
% second gene is the most significant digit, while the last gene 
% is the least significant digit.
%
%   To determine the relative magnitude of a trait, this software uses 
% the DECIMAL[] constant.  This number determines how many digits to the 
% left of the decimal to place the first digit of the trait (=second gene 
% on the trait).
%     Ex:  Suppose trait #i = 123456
%
%     If DECIMAL[i] =       Then Trait[i]=
%
%                 -2                       -0.0023456
%                 -1                       -0.023456
%                  0                       -0.23456
%                  1                       -2.3456
%                  2                      -23.456
%                      etc...
%
% Note also that you must set LOWTRAIT, HIGHTRAIT, DECIMAL, and SIG_FIGS
% so that if we saturate a trait at the value of LOWTRAIT or HIGHTRAIT 
% then the genes on the trait must be able to represent LOWTRAIT or
% HIGHTRAIT.  For example you would not choose HIGHTRAIT(1)=100
% and DECIMAL(1)=-2 with SIG_FIGS(1)=1.  Why?  Hence, you must choose
% the elements of LOWTRAIT and HIGHTRAIT so that they have the same
% decimal position and number of significant figures as their
% corresponding traits.
%
%
% 5.  MUTAT_PROB - the probability that a mutation will occur in any given
%                  gene (i.e., that a gene will be switched to another
%                  element of the alphabet that is being used).  Note that 
%                  for base 10 we randomly choose a number between 0-9, 
%				   with the exception of the number that was being mutated.
%
% 6.  CROSS_PROB - the probability that a crossover will occur between two
%                  parents (standard swapping of genes is used).
%
% 7.  SELF_ENTERED - if this flag is set to 1, an initial 
%                    population can be entered by the user.  If it is 0, 
%                    a random initial population is created.  We enter the 
%                    initial population in a special convenient way for 
%					 the user.
%
% 8. POP_SIZE - the number of individuals in a population (fixed). 
%				 This should be an even number - since pairs mate.
%
% 9. ELITISM - if this flag is set to one, the most fit individual in each
%               generation will be carried over to the next generation
%               without being modified in any way by the genetic operators.
%
% 10. DELTA, EPSILON - these are for the termination criterion.  
%                         The program is terminated if the fitness 
%                         changes less than EPSILON over DELTA
%                         generations.
%
% 11. MAX_GENERATION - the number of generations that the program will 
%                      produce before stopping.  This indicates the length 
%                      of the program will run. If it takes too long 
%                      to run, this number can be lowered.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Next, we provide some ideas on how to run to the program to illustrate
% its operation.  Proceed as follows:
% 
% (1) Experiment with the intial population.  Enter your own initial 
%     population (to do this you modify the code and let SELF_ENTERED = 1).
%     Currently, the code sets the intial population to be all zeros if
%     you let SELF_ENTERED = 1.  Another choice you may want to consider
%     is to place the initial members of the population on some type of
%     uniform grid to try to make sure that at least one of the intial 
%     points is in the proximity of the final point.  Currently, if you let
%     SELF_ENTERED = 0, the program automatically generates a random initial
%     population with a uniform distribution on the allowable range of each 
%     trait.
% (2) Next, you can experiment with the settings for the crossover and
%     mutation probabilities.  First, try MUTAT_PROB=0.05; CROSS_PROB=0.8.
%     If you lower the crossover probability you will generally find that 
%     there will be less local search near the best individuals (i.e., 
%     it will tend to less often interpolate between good individuals to
%     try to find a solution).  If you raise the mutation probability it
%     will tend to cause more random excursions into possibly unexplored 
%     regions.  Too high of a mutation probability can cause degradation
%     to random search, where too low of a mutation probability can cause
%     the GA to lock in on a local minima and never find the global one (at
%     least in a reasonable amount of time).  
% (3) Next, you can experiment with the population size.  Larger population
%     sizes require more computation time but they can speed convergence.
%     For instance, for a larger population size you will have more guesses
%     of where the global maximum is for the intial population and hence
%     it is more likely that there will be at  least one good intial guess.
%     There are similar effects as the algorithm runs.
% (4) Next, turn on elitism.  You will find that elitism removes some of the 
%     effects of having a higher mutation rate since it makes sure that 
%     the best candidate is not disturbed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear           % Initialize momery
rand('state',0) % Reset the random number generator so each time you re-run
                % the program with no changes you will get the same results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User: Input your own parameters below (the ones given are for the
%       optimization problem given above):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  NUM_TRAITS=2;          % Number of traits in each individual
  HIGHTRAIT=[30 30];  	 % Upper limit of a trait
  LOWTRAIT=[0 0]; 	     % Lower limit of a trait
  SIG_FIGS=[6 6]';       % Number of genes in each trait
  DECIMAL=[2 2];         % Order of magnitude the trait
  MUTAT_PROB=0.05;       % Probability of mutation (typically <.1)
  CROSS_PROB=0.8;        % Probability of crossover (typically near 1)
  SELF_ENTERED=0;        % "0": a random initial population.
                         % "1": a specified initial population 
                         % If you choose "1", enter it in the program.
  POP_SIZE=20;           % Number of individuals in the population
  ELITISM=0;             % Elitism ON/OFF, 1/0
  DELTA=100;             % Number of generations to be counted 
                         % for the termination criteria.
  EPSILON = 0.01;  	     % Range that the fitness must change 
                         % in the termination criteria.
  MAX_GENERATION=1000;   % Number of times the loop will run before 
  						 % giving up on the EPSILON-DELTA termination cond.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next, we want to make the initial population of size = POP_SIZE called pop.
% Each column represents a chromosome and each element in that column
% representing a gene.  There will be POP_SIZE chromosomes with each of
% them having CHROM_LENGTH genes.  Now, rather than directly create this
% matrix we allow the user to enter the more easily specified values of the
% traits (pop is coded to be operated on by a genetic algorithm while trait is
% a matrix with each of actual values of the parameters).  Taking this approach 
% you can easily enter an intial population of parameter guesses on the domain
% of the function to be maximized and the first iteration in the GA will
% convert the values to a form (stored in pop) that the genetic operators can 
% work on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

popcount=1;         	% Initialize the generation count, 
						% set it to one, the first population

if SELF_ENTERED == 0            % Make a random initial population for
								% base-10 operation by specifying the 
								% matrix of initial traits
for pop_member = 1:POP_SIZE
      for current_trait = 1:NUM_TRAITS,
 		trait(current_trait,pop_member,popcount)=...
		(rand-(1/2))*(HIGHTRAIT(current_trait)-LOWTRAIT(current_trait))+... 
		(1/2)*(HIGHTRAIT(current_trait)+LOWTRAIT(current_trait));
			    % This starts the population off with numbers chosen randomly 
				% on the allowed range of variation, for this example.  
	  end
end 

else
								
for pop_member = 1:POP_SIZE
        for current_trait = 1:NUM_TRAITS,
            trait(current_trait,pop_member,popcount)=0;  
			% To start with a guess where all the population members are zero				
		end
end 
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next, we need to set some values based on the values input by 
% the user.  In particular, we determine the length of the 
% chromosome and start point of each trait.  This information 
% will be used later in the main algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   CHROM_LENGTH=sum(SIG_FIGS)+NUM_TRAITS; % Length of the chromosome is 
  										  % the number of sig. figs. plus
  										  % the number of sign positions
   TRAIT_START(1)=1;					  % Initialize: the first trait 
   										  % starts at the first digit 
   										  % (this is the sign digit)

   for current_trait=1:NUM_TRAITS, 	% Determine the start point of the
   									% other traits - it is the start of
   									% the last trait plus the no. of sig. 
   									% figs. plus one for sign
TRAIT_START(current_trait+1)=...
   TRAIT_START(current_trait)+SIG_FIGS(current_trait)+1;
   % Yes, we compute the TRAIT_START for one extra trait - this 
   % is used for convenience in the code below.
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this next loop the fitness is calculated, the children are 
% created and it repeats until the EPSILON-DELTA termination condition
% is satisfied or MAX_GENERATION is reached.  This is the main
% loop of the GA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while popcount <= MAX_GENERATION 
 
% First, fix bad traits (i.e., ones that are out of the range 
% specified by HIGHTRAIT and LOWTRAIT) by saturation at the extremes

for pop_member = 1:POP_SIZE

for current_trait = 1:NUM_TRAITS,


if trait(current_trait,pop_member,popcount)>HIGHTRAIT(current_trait)

                     % The trait has went higher than the upper
                     % bound so let the trait equal to the
                     % HIGHTRAIT bound.

trait(current_trait,pop_member,popcount)=HIGHTRAIT(current_trait);

% Now consider the other case:

elseif trait(current_trait,pop_member,popcount)<LOWTRAIT(current_trait)

                     % The trait has went lower than the lower 
                     % bound so let the trait equal to the
                     % LOWTRAIT bound

trait(current_trait,pop_member,popcount)=LOWTRAIT(current_trait);
                  
		end
               
% Now that we have reset the traits to be in range, we must 
% convert them to the chromosome form for use with the genetic operators.
% First, we transfer the sign of the trait into the chromosome

               if trait(current_trait,pop_member,popcount) < 0
                  pop(TRAIT_START(current_trait),pop_member)=0;
               else
                  pop(TRAIT_START(current_trait),pop_member)=9;
               end

% Next, strip off the sign and store the resulting value in a
% temporary variable that is used in the construction of pop

temp_trait(current_trait,pop_member)=...
           abs(trait(current_trait,pop_member,popcount));  
		% temp_trait is trait without the sign of trait 

% Next, we store the numbers of the trait in the chromosome:
% First, set up a temporary trait with at most
% one nonzero digit to the left of the decimal point.
% This is used to strip off the numbers to put 
% them into a chromosome.

temp_trait(current_trait,pop_member)=...
temp_trait(current_trait,pop_member)/10^(DECIMAL(current_trait)-1);

% Encode the new trait into chromosome form 

for make_gene = TRAIT_START(current_trait)+1:TRAIT_START(current_trait+1)-1,
               
% For each gene on the trait make the gene the corresponding digit on 
% temp_trait (note that rem(x,y)=x-roundtowardszero(x/y)*y or 
% rem(x,1)=x-roundtowardszero(x) so that rem(x,1) is the fraction part
% and x-rem(x,1) is the integer part so the next line makes the location in
% pop the same as the digit to the left of the decimal point of temp_trait
  
pop(make_gene,pop_member)=temp_trait(current_trait,pop_member)-...
   rem(temp_trait(current_trait,pop_member),1);

% Next, we take temp_trait and rotate the next digit to the left so that
% next time around the loop it will pull that digit into the 
% chromosome.  To do this we strip off the leading digit then shift 
% in the next one.
  
temp_trait(current_trait,pop_member)=...
  (temp_trait(current_trait,pop_member)-pop(make_gene,pop_member))*10;
               
end

         end % Ends "for current_trait=..." loop
         
	 end    % Ends "for pop_member=..." loop
	 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER: Below is where you change the fitness function that is to be
% maximized.  This may involve changing several lines of code below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the example below we want to *maximize* the function z=f(x,y).
% Important: When defining the fitness function you must remember that 
% it is a function that you are *maximizing* and that it must always be 
% positive (the selection mechanism depends on this).  For our example, 
% we are maximizing the function f.  Below, we include some ideas about
% what to do if you seek to minimize a function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sumfitness = 0; 		% Re-initialize for each generation

% First, determine the values of the function to be minimized for 
% each chromosome.

   for chrom_number = 1:POP_SIZE,       % Test fitness

fitness_bar(chrom_number)=...
+5*exp(-0.1*((trait(1,chrom_number,popcount)-15)^2+...
   (trait(2,chrom_number,popcount)-20)^2))...
-2*exp(-0.08*((trait(1,chrom_number,popcount)-20)^2+...
   (trait(2,chrom_number,popcount)-15)^2))...
+3*exp(-0.08*((trait(1,chrom_number,popcount)-25)^2+...
   (trait(2,chrom_number,popcount)-10)^2))...
+2*exp(-0.1*((trait(1,chrom_number,popcount)-10)^2+...
   (trait(2,chrom_number,popcount)-10)^2))...
-2*exp(-0.5*((trait(1,chrom_number,popcount)-5)^2+...
   (trait(2,chrom_number,popcount)-10)^2))...
-4*exp(-0.1*((trait(1,chrom_number,popcount)-15)^2+...
   (trait(2,chrom_number,popcount)-5)^2))...
-2*exp(-0.5*((trait(1,chrom_number,popcount)-8)^2+...
   (trait(2,chrom_number,popcount)-25)^2))...
-2*exp(-0.5*((trait(1,chrom_number,popcount)-21)^2+...
   (trait(2,chrom_number,popcount)-25)^2))...
+2*exp(-0.5*((trait(1,chrom_number,popcount)-25)^2+...
   (trait(2,chrom_number,popcount)-16)^2))...
+2*exp(-0.5*((trait(1,chrom_number,popcount)-5)^2+...
   (trait(2,chrom_number,popcount)-14)^2));

   end

% Next, compute the fitness function (this loop is only kept 
% separate from the next one in case you need to implement 
% some of the options discussed in the comments).  

   for chrom_number = 1:POP_SIZE,       % Test fitness

% The fitness function must be chosen so that it is always positive.
% To ensure this for our example we assume that we know that the value
% of f will never be below -5; so we simply add 5 to every fitness_bar value.
% Another approach would be to simply find the minimum value 
% of the fitness_bar function for every element in the population 
% and then add on its absolute value so that all the resulting 
% values will be positive

 fitness(chrom_number)= fitness_bar(chrom_number) + 5;
 
% Notice that when we turn a minimization problem into a maximization
% problem we often use 
% fitness(chrom_number)= -fitness_bar(chrom_number) + max(fitness_bar);
% as the above line.  Here, the minus sign is
% for turning the minimization problem into a maximization problem
% and max(fitness_bar) is used to make the fitness function positive.
% Another option would be to use the fitness function:
% fitness(chrom_number)=(1/(fitness_bar(chrom_number) + .1));
% The 1/fitness_bar function switches it to a maximzation problem.  
% The +.1 is to make sure that there is no divide by zero. Note that you can
% tune the GA in this case by changing the 1 in the numerator to another 
% constant and the .1 to a different value.

sumfitness = sumfitness + fitness(chrom_number);  % Store this for 
      											  % use below
   end

% Next, determine the most fit and least fit chromosome and 
% the chrom_numbers (which we call bestmember and worstmember). 

[bestfitness(popcount),bestmember]=max(fitness);
[worstfitness(popcount),worstmember]=min(fitness);

% Next, save these (if want to save worstindividual can too)

bestindividual(:,popcount)=trait(:,bestmember,popcount);
%worstindividual(:,popcount)=trait(:,worstmember,popcount);

% Compute the average fitness in case you want to plot it.

avefitness(popcount) = sumfitness / POP_SIZE;        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the next generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, form the mating pool.  
% To do this we select as parents the
% chromosomes that are most fit.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for pop_member = 1:POP_SIZE,
	if ELITISM ==1 & pop_member==bestmember  % If elitism on, and have
											 % the elite member
		parent_chrom(:,pop_member)=pop(:,pop_member); % Makes sure that
		                          % the elite member gets into the next
								  % generation.
	else
      pointer=rand*sumfitness; % This makes the pointer for the roulette 
                               % wheel.  
      member_count=1;        % Initialization
      total=fitness(1);

      while total < pointer,                % This spins the wheel to the 
       									    % pointer and finds the 
       									    % chromosome there - which is
       									    % identified by member_count
         member_count=member_count+1;
         total=total+fitness(member_count);
      end
      
% Next, make the parent chromosome

    parent_chrom(:,pop_member)=pop(:,member_count);
	end
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproduce section (i.e., make off-spring - "children")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the approach below each individual gets to mate and
% they randomly pick someone else (not themselves) in the 
% mating pool to mate with.  Resulting children are 
% composed of a combination of genetic material of their
% parents.  If elitism is on, when the elite member gets
% a chance to mate they do not; they are simply copied
% over to the next generation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for parent_number1 = 1:POP_SIZE,    % Crossover (parent_number1 is the 
									% individual who gets to mate
	if ELITISM ==1 & parent_number1==bestmember % If elitism on, and 
											 % have the elite member
		child(:,parent_number1)=parent_chrom(:,parent_number1);
	else
	   parent_number2=parent_number1;  % Initialize who the mate is
	   while parent_number2 == parent_number1   % Iterate until find 
		                                        % a mate other than
												% yourself
         parent_number2 = rand*POP_SIZE; % Choose parent number 2
		 								 % randomly (a random mate)
         parent_number2 = parent_number2-rem(parent_number2,1)+1;     
	   end
    if CROSS_PROB > rand           % If true then crossover occurs

         site = rand*CHROM_LENGTH;   % Choose site for crossover
         site = site-rem(site,1)+1;  % and make it a valid integer 
         							 % number for a site

% The next two lines form the child by the swapping of genetic
% material between the parents

child(1:site,parent_number1)=parent_chrom(1:site,parent_number1);

child(site+1:CHROM_LENGTH,parent_number1)=...
     parent_chrom(site+1:CHROM_LENGTH,parent_number2);

      else                           % No crossover occurs

% Copy non-crossovered chromosomes into next generation
% In this case we simply take one parent and make them
% the child.

child(:,parent_number1)=parent_chrom(:,parent_number1);

	 end
   end  % End the "if ELITISM..." statement
end  % End "for parent_number1=..." loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mutate children.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, we mutate to a different allele
% with a probability MUTAT_PROB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for pop_member= 1:POP_SIZE,

	if ELITISM ==1 & pop_member==bestmember  % If elitism on, and 
											 % have the elite member
		child(:,pop_member)=child(:,pop_member); % Do not mutate
											% the elite member
	else
	  for site = 1:CHROM_LENGTH,

         if MUTAT_PROB > rand        % If true then mutate  
            rand_gene=rand*10; 		 % Creat a random gene
            
            % If it is the same as the one already there then 
            % generate another random allele in the alphabet
            
            while child(site,pop_member) == rand_gene-rem(rand_gene,1),
               rand_gene=rand*10;
            end;
            
            % If it is not the same one, then mutate
            
            child(site,pop_member)=rand_gene-rem(rand_gene,1);
            
            % If takes a value of 10 (which it cannot
            % mutate to) then try again (this is a very low probability
			% event (most random number generators generate numbers
			% on the *closed* interval [0,1] and this is why this line
			% is included).
            
            if rand_gene == 10
                 site=site-1;
            end
         end  % End "if MUTAT_PROB > rand ... 
	   end  % End for site... loop
	 end  % End "if ELITISM..."
   end  % End for pop_member loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the next generation (this completes the main part of the GA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   pop=child;               % Create next generation (children 
   							% become parents)       

   popcount=popcount+1;		% Increment to the next generation
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next, we have to convert the population (pop) to the base-10 
% representation (called trait) so that we can check if the traits
% all still lie in the proper ranges specified by HIGHTRAIT and LOWTRAIT
% at the beginning of the next time around the loop.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pop_member = 1:POP_SIZE

	for current_trait = 1:NUM_TRAITS,
			 
trait(current_trait,pop_member,popcount)=0; % Initialize variables
place_pointer=1;

% Change each of the coded traits on the chromosomes into base-10 traits: 
% For each gene on the current_trait past the sign digit but before the
% next trait find its real number amount and hence after finishing
% the next loop trait(current_trait,pop_member,popcount) will be the base-10
% number representing the trait

for gene=TRAIT_START(current_trait)+1:TRAIT_START(current_trait+1)-1,
 place=DECIMAL(current_trait)-place_pointer;
 trait(current_trait,pop_member,popcount)=...
 trait(current_trait,pop_member,popcount)+...
    (pop(gene,pop_member))*10^place;
 place_pointer=place_pointer+1;
end

% Determine sign of the traits and fix 
% trait(current_trait,pop_member,popcount) so that it has the right sign:

if pop(TRAIT_START(current_trait),pop_member) < 5
   trait(current_trait,pop_member,popcount)=...
      -trait(current_trait,pop_member,popcount);
end


	end % Ends "for current_trait=..." loop
         
end    % Ends "for pop_member=..." loop
  			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terminate the program when the best fitness has not changed
% more than EPSILON over the last DELTA generations.  It would also
% make sense to use avefitness rather than bestfitness in this test.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if popcount > DELTA+1 & ...
       max(abs(bestfitness(popcount-DELTA:popcount-1)-...
          bestfitness(popcount-DELTA-1:popcount-2)))<=EPSILON
       break;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end  % End "for pop_count=..." loop - the main loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next, provide some plots of the results of the simulation.
% These are for the example function that we seek to maximize; however
% for other applications it is not difficult to modify this code
% to show plots that are informative for the application at hand.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=1:popcount-1;  % For use in plotting

x=0:31/100:30;   % For our function the range of values we are considering
y=x;

% Compute the function that we are trying to find the maximum of.

for jj=1:length(x)
	for ii=1:length(y)
		z(ii,jj)=...
		+5*exp(-0.1*((x(jj)-15)^2+(y(ii)-20)^2))...
		-2*exp(-0.08*((x(jj)-20)^2+(y(ii)-15)^2))...
		+3*exp(-0.08*((x(jj)-25)^2+(y(ii)-10)^2))...
		+2*exp(-0.1*((x(jj)-10)^2+(y(ii)-10)^2))...
		-2*exp(-0.5*((x(jj)-5)^2+(y(ii)-10)^2))...
		-4*exp(-0.1*((x(jj)-15)^2+(y(ii)-5)^2))...
		-2*exp(-0.5*((x(jj)-8)^2+(y(ii)-25)^2))...
		-2*exp(-0.5*((x(jj)-21)^2+(y(ii)-25)^2))...
		+2*exp(-0.5*((x(jj)-25)^2+(y(ii)-16)^2))...
		+2*exp(-0.5*((x(jj)-5)^2+(y(ii)-14)^2));
	end
end

% First, show the actual function to be maximized and its contour map

figure(1)
clf
surf(x,y,z);
colormap(jet)
% Use next line for generating plots to put in black and white documents.
colormap(white);
xlabel('x');
ylabel('y');
zlabel('z');
title('Fitness function');

figure(2) 
clf
contour(x,y,z,25)
colormap(jet)
% Use next line for generating plots to put in black and white documents.
colormap(gray);
xlabel('x');
ylabel('y');
zlabel('z');
title('Fitness function (contour map)');

hold on

% Next, add on the evolution of the other population
% members over time.

for pop_member = 1:POP_SIZE
	xx=squeeze(trait(1,pop_member,:));
	yy=squeeze(trait(2,pop_member,:));
	traitplot=plot(xx,yy,'k.');
	set(traitplot,'MarkerSize',6);
end

% Then, add on the evolution of the best individual over 
% time that the GA finds.

btraitplot=plot(bestindividual(1,:),bestindividual(2,:),'ko');
set(btraitplot,'MarkerSize',4)

hold off

% Next, plot some data on the operation of the GA

figure (3)
clf
subplot(211)
plot(t,bestfitness,'k--',t,avefitness,'k-',t,worstfitness,'k-.')
xlabel('Generation')
ylabel('Best, average, and worst fitness')
title('Best and worst fitness vs. generation')
subplot(212)
plot(t,bestindividual(1,:),'k-',t,bestindividual(2,:),'k--')
xlabel('Generation')
ylabel('Best individuals')
title('Best individuals vs. generation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
