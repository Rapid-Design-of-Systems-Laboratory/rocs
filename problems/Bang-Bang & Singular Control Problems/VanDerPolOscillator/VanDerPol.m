function [in] = VanDerPol() 
% 
% This function creates input for the optimization problem,  
% which is later used by other functions  
% 
% input : void 
% output : in [structure] 
% Developed by : Dr. K. Mall
% Last modified: Jan 29, 2020

%%%%%%%%%%%%%%%%%%%%%%% 
%% Execution Control %% 
%%%%%%%%%%%%%%%%%%%%%%% 

% Set desired tolerances for the solvers 
in.oc.writeEquations = true;
in.useDeval          = true;

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%

% Doubles must be equal to one (not dynamically updated during continuation)
in.autoScale = true;
in.scale = {'m',1; ...
			's',1; ...
			'nd',1}; % nd = nondimensional

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Independent Variable %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g., time

in.oc.independentVariable = {'t','s'}; % time

%%%%%%%%%%%%
%% States %%
%%%%%%%%%%%%

in.oc.state = {'x1','nd'; ... % radial position magnitude
			   'x2','nd';...
               'x3','nd';...
               't','s'}; % relative velocity

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%

in.oc.stateRate = {'x2'; ...
				   '-x1 + x2*(1-x1^2) + sin(vdpcont)';...
                   '0.5*(x1^2 + x2^2)';...
                   '1'}; 

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control = {'vdpcont','nd'};

in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path cost
in.oc.cost.path     = {'epsilon*cos(vdpcont)','nd'}; 

% Terminal cost 
in.oc.cost.terminal = {'x3','nd'}; 

% Initial cost 
in.oc.cost.initial  = {'0','nd'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial constraint 
in.oc.constraint.initial = {'x1-x0(1)','nd'; ...
							'x2-x0(2)','nd';...
                            'x3-x0(3)','nd';...
                            't-x0(4)','s'};
													 
% Terminal constraint
in.oc.constraint.terminal = {'t-xf(4)','s'};

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%

% in.oc.initialGuessFunc = @getInitialGuess;
in.oc.guess.mode = 'auto';
in.oc.guess.timeIntegrate = 1.5; % 0.1 leads to local min

% % Use automatic init
% Conditions at entry
in.oc.guess.initial.x1 = 0;
in.oc.guess.initial.x2 = 1;
in.oc.guess.initial.x3 = 0;
in.oc.guess.initial.t = 0;
in.oc.costates = -0.1;

% Constants
in.const.epsilon  = {0.1,'nd'};% {0.1,'nd'}; % scaling factor for smoothed bang-bang control
in.const.tol      = {1e-6,'nd'}; 
in.const.NMax     = {1e13,'nd'}; 
in.const.NMesh    = {500,'nd'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Execution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = 0;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%%

ind = ind+1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 10;% Number of steps in the continuation set 

in.CONT{ind}.constraint.terminal.t = 4;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%%run
ind = ind + 1;
in.cont.method(ind) = 1;

in.CONT{ind}.numCases = 200;% Number of steps in the continuation set
in.CONT{ind}.const.epsilon = linspace(0,-0.0999,in.CONT{ind}.numCases);

return