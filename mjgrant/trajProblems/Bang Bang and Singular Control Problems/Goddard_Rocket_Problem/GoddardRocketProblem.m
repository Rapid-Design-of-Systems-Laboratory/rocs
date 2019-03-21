function [in] = GoddardRocketProblem()
%
% This function creates input for the optimization problem,
% which is later used by other functions
%
% input : void
% output : in [structure]
% Developed by : Dr. M.J. Grant, Kshitij Mall and Thomas Antony
% Last modified: 25 Mar, 2014

%%%%%%%%%%%%%%%%%%%%%%%
%% Execution Control %%
%%%%%%%%%%%%%%%%%%%%%%%
in.oc.writeEquations = false; % Determine if we need to regenerate the equation files
in.bvpOptions = bvpset('AbsTol',1e-6,'RelTol',1e-6,'Nmax',100000000,'Stats','on'); 

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%
% These scaling parameters are altered automatically during continuation. This
% is just a first guess.
% Doubles must be equal to one (not dynamically updated during continuation)
in.autoScale = true;

in.scale = {'ft',1; ...
		  'rad',1; ...
		  's',1; ...
		  'lbm',1; ...
		  'nd',1}; % nd = nondimensional

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Independent Variable %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g., time

in.oc.independentVariable = {'t','s'}; % time

%%%%%%%%%%%%
%% States %%
%%%%%%%%%%%%

% States
in.oc.state = {'h','ft'; ... % altitude
		 'v','ft/s'; ... % relative velocity
		 'mass','lbm'}; %... % relative flight-path angle ...
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%
% T = '(96.5*sin(control)+96.5)';
T = '(0.5*TMax*(sin(control)+1))';

% Quantities of Interest
D = 'dragk*v^2*exp(-h/H)'; % Drag Force [N]

% Equations of Motion
%% Planar
in.oc.stateRate = {'v + epsilon*cos(control)'; ...
				['(',T,'-',D,')/mass - g0']; ...
				['-',T,'/c']};

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control = {'control','nd'}; % angle of attack control
in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximize/minimize
in.minimize = true;

% Path cost
in.oc.cost.path = {'0','ft'};

% Terminal cost
in.oc.cost.terminal = {'-h','ft'};

% Initial cost
in.oc.cost.initial = {'0','ft'};

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial constraint
in.oc.constraint.initial = {'h-x0(1)','ft'; ...
	 						 'v-x0(2)','ft/s'; ...
	 						 'mass-x0(3)','lbm'};

% Terminal constraint
in.oc.constraint.terminal = {'mass-xf(3)','lbm'};

%%%%%%%%%%%%%%%%%%%%%%
%%     Constants    %%
%%%%%%%%%%%%%%%%%%%%%%
in.const.dragk = {5.49153485e-5,'nd'}; % Gravitational parameter, m^3/s^2
in.const.rho0 = {0.002378,'lbm/ft^3'}; % Sea-level atmospheric density, kg/m^3
in.const.H = {23800,'ft'}; % Scale height for atmosphere of Earth, m
in.const.c = {1.580942527987656e+03,'ft/s'}; % Radius of planet, m
in.const.g0 = {32.174,'ft/s^2'};
in.const.epsilon = {1,'s/lbm'};
in.const.TMax = {193,'lbm*ft/s^2'};

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%
% in.oc.initialGuessFunc = @getInitialGuess;

in.oc.guess.mode = 'auto';
in.oc.guess.timeIntegrate = 20; % 0.1 leads to local min

% % Use automatic init
% Conditions at entry
in.oc.guess.initial.h = 0;
in.oc.guess.initial.v = 0;
in.oc.guess.initial.mass = 3;
in.oc.costates = -0.1;

%%%%%%%%%%%%%%%%%%
%% Continuation %%
%%%%%%%%%%%%%%%%%%

in.cont.method = 1; % 1 = manually changing parameters

ind = 0;
%
ind = ind+1;
in.CONT{ind}.numCases = 10; % 41; Number of steps in the continuation set
in.CONT{ind}.constraint.terminal.mass = 1;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%%
ind = ind + 1;
in.cont.method(ind) = 1;

in.CONT{ind}.numCases = 10; %96;%9; % Number of steps in the continuation set
% in.CONT{ind}.const.epsilon = 0.95*linspace(0,-0.999,in.CONT{ind}.numCases);
in.CONT{ind}.const.epsilon = linspace(0,-0.5,in.CONT{ind}.numCases);

return
