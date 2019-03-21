function [in] = defaultConfig() 
% 
% This function creates default configuration options for the optimization problem,  
% which is later used by other functions  
% 
% input : void 
% output : in [structure] 
% Developed by : Thomas Antony
% Created : 13 January, 2015

%%%%%%%%%%%%%%%%%%%%%%% 
%% Execution Control %% 
%%%%%%%%%%%%%%%%%%%%%%% 
 
in.run.optimalCalcs = true; % Calculate necessary conditions of optimality 
in.run.continuation = true; % Execute continuation method on the shorter solution 
in.autocodeDirTraj = [pwd,'/autocode']; % location of autocoded functions 

% When verbose is enabled on warnings, MATLAB displays an extra line of  
% information with each warning that tells about way to supress it  

in.verbose = false;
in.writeJac = false; % Write Jacobian information in BVP4C 
in.useJac = false; % Use Jacobian information in BVP4C 
in.useMex = true;
in.skipUnconvergedSolutions = false;

% Continuation variables (ignore the 'in' part since common for all variables) 
in.continuationVarSet = {'const','constraintVal'}; 
in.convertParametersToStates = false; % Consider parameters like tf and 'nu's as extra "states" 

in.rootSolving = 0; % 0 = bvp4c, 1 = GPU STM shooting 
% Set desired tolerances for the solvers 
in.odeOptions = odeset('AbsTol',1e-8,'RelTol',1e-8); 
in.bvpOptions = bvpset('AbsTol',1e-4,'RelTol',1e-4,'Nmax',100000,'Stats','on'); 

in.gpuSolve.eomFormat = 'cubin'; 
in.gpuSolve.sttOrder  = 1; 
in.gpuSolve.TimeSteps = 1024; 	% This number must be an integer multiple of 

in.gpuSolve.ArcSegments = 256; 
in.gpuSolve.Tol = 1e-4; 
in.gpuSolve.dampingFactor = 1.00;
in.gpuSolve.maxIterations = 200;
in.gpuSolve.write = true;

in.maxNumArcs = 1; % maximum number of arcs for variable size array allocation during mex-ing

% Determine if write STT function. Only need to do this when order of STT 
% changes. 
in.oc.writeSTTfunc = false; 
in.oc.writeEquations = true;

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%
% These scaling parameters are altered automatically during continuation. This 
% is just a first guess. 
% Doubles must be equal to one (not dynamically updated during continuation)
in.autoScale = true;

in.scale = {'m',1; ...
                      'rad',1; ...
                      's',1; ...
                      'kg',1; ...
                      'nd',1}; % nd = nondimensional
 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Independent Variable %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g., time
in.oc.independentVariable = {'t','s'}; % time
            
%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%
in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximize/minimize
in.minimize = true;

% Path cost
in.oc.cost.path = {'0','m^2/s^2'}; 

% Terminal cost 
in.oc.cost.terminal = {'0','m^2/s^2'}; 

% Initial cost 
in.oc.cost.initial = {'0','m^2/s^2'}; 

in.oc.cost.interiorPoint = {'','',0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In-Flight Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% 
%%     Constants    %% 
%%%%%%%%%%%%%%%%%%%%%% 

% Smoothing parameter for bang-bang problems
in.const.epsilon = {1e-2,'nd'};

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%
in.oc.guess.timeIntegrate = 0.1;
in.oc.guess.costate = -0.1;
% in.oc.initialGuessFunc = @getInitialGuess;

%%%%%%%%%%%%%%%%%% 
%% Continuation %% 
%%%%%%%%%%%%%%%%%% 

in.cont.method = 1; % 1 = manually changing parameters, 2 = targeting new solution using STTs 
in.cont.orderSTT = 0; % Order of state transition tensors used to estimate new solutions 

in.cont.indexContBeforeTrades = inf; % This ensures constraint values not updated after for future continuations (for trades in constraints themselves, etc.) May be able to delete this in the future.

in.cont.startIndex = 1;

return

