function [in] = secondorder() 
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

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%

% Doubles must be equal to one (not dynamically updated during continuation)
in.autoScale = true;
in.scale = {'m',1; ...
			's',1; ...
            'rad',1;...
			'nd',1}; % nd = nondimensional

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Independent Variable %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g., time

in.oc.independentVariable = {'t','s'}; % time

%%%%%%%%%%%%
%% States %%
%%%%%%%%%%%%

in.oc.state = {'x','m'; ... % position 
			   'v','m/s';...% velocity
               'z','s'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%

in.oc.stateRate = {'v'; ...
				   'umax*sin(socontrol)';...
                   '1'}; 

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control = {'socontrol','rad'};

in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path cost
in.oc.cost.path = {'umax^2*sin(socontrol)^2/2','m^2/s^4'}; 

% Terminal cost 
in.oc.cost.terminal = {'0','m^2/s^4'}; 

% Initial cost 
in.oc.cost.initial = {'0','m^2/s^4'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial constraint 
in.oc.constraint.initial = {'x-x0(1)','m'; ...
							'v-x0(2)','m/s';...
                            'z-x0(3)','s'};
													 
% Terminal constraint
in.oc.constraint.terminal = {'x-xf(1)','m'; ...
							'v-xf(2)','m/s';...

                            'z-xf(3)','s'};
% Constants
in.const.umax    = {3,'m/s^2'};    
in.const.tol     = {1e-4,'nd'}; 
in.const.NMax    = {1e10,'nd'}; 

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%

% in.oc.initialGuessFunc = @getInitialGuess;
in.oc.guess.mode = 'auto';
in.oc.guess.timeIntegrate = 2; % 0.1 leads to local min

% % Use automatic init
% Conditions at entry
in.oc.guess.initial.x = 1;
in.oc.guess.initial.v = 1;
in.oc.guess.initial.z = 0;
in.oc.costates = -0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Execution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = 0;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%%

ind = ind+1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 10;% 7 Number of steps in the continuation set 

in.CONT{ind}.constraint.terminal.x = 0;
in.CONT{ind}.constraint.terminal.v = 0;
in.CONT{ind}.constraint.terminal.z = 2;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%%
ind = ind + 1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 10;% Number of steps in the continuation set 
in.CONT{ind}.const.umax = linspace(0,-0.5,in.CONT{ind}.numCases);

return