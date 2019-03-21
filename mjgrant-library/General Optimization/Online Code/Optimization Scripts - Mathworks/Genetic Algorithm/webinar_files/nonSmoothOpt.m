%% Genetic Algorithm optimization demo
% This is a demonstration of how to find a minimum of a non-smooth
% objective function using the Genetic Algorithm (GA) function in the
% Genetic Algorithm and Direct Search Toolbox. Traditional derivative-based
% optimization methods, like those found in the Optimization Toolbox, are
% fast and accurate for many types of optimization problems.  These methods
% are designed to solve 'smooth', i.e., continuous and differentiable,
% minimization problems, as they use derivatives to determine the direction
% of descent. While using derivatives makes these methods fast and
% accurate, they often are not effective when problems lack smoothness,
% e.g., problems with discontinuous, non-differentiable, or stochastic
% objective functions. When faced with solving such non-smooth problems,
% methods like the genetic algorithm or the more recently developed pattern
% search methods, both found in the Genetic Algorithm and Direct Search
% Toolbox, are effective alternatives. 

%% Initialization
clear all; close all;format compact
Objfcn = @nonSmoothFcn;   %Handle to the objective function
X0 = [2 -2];   % Starting point 
range = [-6 6;-6 6];      %Range used to plot the objective function
rand('state',0);          %Reset the state of random number generators
randn('state',0);

%% Non-smooth Objective Function
% To help visualize the problem and results, we have chosen a problem with
% only two variables, but the algorithms we explore are certainly not
% limited to such small problem. We can view the code for this objective
% function.
type nonSmoothFcn.m
%%
% The objective function of our sample problem is a piece-wise continuous
% function, that is it has smooth regions separated by discontinuities,
% with one exceptional region that is non-differentiable almost everywhere.
% We use the function PLOTOBJ, written for this demo, to plot the function
% NONSMOOTHFCN over the range = [-6 6;-6 6].
showNonSmoothFcn(Objfcn,range);
set(gca,'CameraPosition',[-36.9991   62.6267  207.3622]);
set(gca,'CameraTarget',[0.1059   -1.8145   22.3668])
set(gca,'CameraViewAngle',6.0924)
%Plot of the starting point (used by the PATTERNSEARCH solver)
plot3(X0(1),X0(2),feval(Objfcn,X0),'or','MarkerSize',10,'MarkerFaceColor','r');
fig = gcf;

%% Minimization Using The Genetic Algorithm
% The motivation for the genetic algorithm is evolutionary biology and
% genetics, mainly Darwin’s theory of survival of the fittest. The Genetic
% Algorithm (GA) works on a population using a set of operators that are
% applied to the population. A population is a set of points in the design
% space. The initial population is generated randomly by default. The next
% generation of the population is computed using the fitness of the
% individuals in the current generation. The genetic algorithm does not use
% derivatives to detect descent in its minimization steps and so it is a
% good choice for problems such as our non-differentiable problem, as well
% as discontinuous, and stochastic problems. 
%%
% To start, we will use the Genetic Algorithm, GA, alone to find the
% minimum of the objective function (also called a fitness function). We
% need to supply GA with a function handle to the fitness function
% nonSmoothFcn.m. Also, GA needs to know the how many variables are in the 
% problem, which is two for this function.
FitnessFcn = @nonSmoothFcn;
numberOfVariables = 2;
%%
% Run GA Alone
% Some plot functions are selected to monitor the performance of the solver.
optionsGA = gaoptimset('PlotFcns',@gaplotbestfun,'PlotInterval',5, ...
                   'PopInitRange',[-5;5]);
% We run GA with the options 'optionsGA' as the third argument.
[Xga,Fga] = ga(FitnessFcn,numberOfVariables,optionsGA)
% Plot the final solution
figure(fig)
hold on;
plot3(Xga(1),Xga(2),Fga,'vm','MarkerSize',10,'MarkerFaceColor','m');
hold off;
fig = gcf;

%%
% The optimum is at x* = (-4.7124, 0.0). GA found the point
% (-4.7775,0.0481) near the optimum, but could not get closer with the
% default stopping criteria. By changing the stopping criteria, we might
% find a more accurate solution, but it may take many more function
% evaluations to reach x* = (-4.7124, 0.0). Instead, we can use a more
% efficient local search that starts where GA left off. The hybrid function
% field in GA provides this feature automatically.

%% Minimization Using A Hybrid Function
% We will use a hybrid function to solve the optimization problem, i.e.,
% when GA stops (or you ask it to stop) this hybrid function will start
% from the final point returned by GA. Our choices are FMINSEARCH,
% PATTERNSEARCH, or FMINUNC. Since this optimization example is smooth near
% the optimizer, we can use the FMINUNC function from the Optimization
% toolbox as our hybrid function as this is likely to be the most
% efficient. Since FMINUNC has its own options structure, we provide it as
% an additional argument when specifying the hybrid function.  

%% 
% Run GA-FMINUNC Hybrid
optHybrid = gaoptimset(optionsGA,'Generations',15, 'PlotInterval',1,...
                     'HybridFcn',{@fminunc,optimset('OutputFcn',@fminuncOut)});
[Xhybrid,Fhybrid] = ga(Objfcn,2,optHybrid);
% Plot the final solution
figure(fig);
hold on;
plot3(Xhybrid(1),Xhybrid(2),Fhybrid+1,'^c','MarkerSize',10,'MarkerFaceColor','c');
hold off;
%%
% The plot shows the best value of the population in every generation. The
% best value found by GA when it terminated is also shown in the plot
% title. When GA terminated, FMINUNC (the hybrid function) was
% automatically called with the best point found by GA so far. The solution
% 'Xhybrid' with function value 'Fhybrid' is the result of using GA and
% FMINUNC together. As shown here, using the hybrid function can
% efficiently improve the accuracy of the solution. Also, the total number
% of generations was reduced from 100 to 18 by using FMINUNC solver as
% hybrid a function. The improvemnet in X and function value is calculated
% below. 
disp(['The norm of |Xga - Xhb| is  ', num2str(norm(Xga-Xhybrid))]);
disp(['The difference in function values Fga and Fhb is ', num2str(Fga - Fhybrid)]);

%% Minimization Using The Pattern Search Algorithm
% Pattern search, although much less well known, is an attractive
% alternative to the genetic algorithm as it is often computationally less
% expensive and can minimize the same types of functions. Additionally, the
% toolbox includes a pattern search method that can solve problems with
% linear constraints and bounds.

%% 
% Pattern search operates by searching a set of points called a pattern,
% which expands or shrinks depending on whether any point within the
% pattern has a lower objective function value than the current point.  The
% search stops after a minimum pattern size is reached. Like the genetic
% algorithm, the pattern search algorithm does not use derivatives to
% determine descent, and so works well on non-differentiable, stochastic,
% and discontinuous objective functions.  And similar to the genetic
% algorithm, pattern search is often very effective at finding a global
% minimum because of the nature of its search. 

%%
% To minimize our objective function using the PATTERNSEARCH function, 
% we need to pass in a function handle to the objective function as well 
% as specifying a start point as the second argument. 
ObjectiveFunction = @nonSmoothFcn;
X0 = [2 -2];   % Starting point 
% Some plot functions are selected to monitor the performance of the solver.
optionsPS = psoptimset('PlotFcns',@psplotbestf);
% Run pattern search solver
[Xps,Fps] = patternsearch(Objfcn,X0,[],[],[],[],[],[],optionsPS)
% Plot the final solution
figure(fig)
hold on;
plot3(Xps(1),Xps(2),Fps+1,'*y','MarkerSize',14,'MarkerFaceColor','y');
hold off;
