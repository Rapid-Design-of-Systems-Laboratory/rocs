function [in] = Aerocapture_bounded()
%
% This function creates input for the optimization problem,
% which is later used by other functions
%
% input : void
% output : in [structure]
% Developed by : Dr. K. Mall, Dr. M.J. Grant, and Dr. T. Antony
% Last modified: Jan 28, 2020

%%%%%%%%%%%%%%%%%%%%%%%
%% Execution Control %%
%%%%%%%%%%%%%%%%%%%%%%%

in.oc.writeEquations = true; % Determine if we need to regenerate the equation files
in.useDeval          = false;

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%
% These scaling parameters are altered automatically during continuation. This
% is just a first guess.
% Doubles must be equal to one (not dynamically updated during continuation)

in.autoScale = true;
in.scale = {'m','x1'; ... % x1 is the first state, altitude
		    'rad',1; ...
		    's','x1/x3'; ... % x3 is the third state, velocity. x1/x3 = time
		    'kg','const.mass'; ...
		    'nd',1}; % nd = nondimensional
      
% in.scale = {'m',1; ...
% 		    'rad',1; ...
% 		    's',1; ...
% 		    'kg',1; ...
% 		    'nd',1}; % nd = nondimensional

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Independent Variable %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g., time

in.oc.independentVariable = {'t','s'}; % time

%%%%%%%%%%%%
%% States %%
%%%%%%%%%%%%

% States
in.oc.state = {'h','m'; ... % altitude
		       'thetta','rad'; ... % longitude, positive eastward
		       'v','m/s'; ... % relative velocity
		       'gam','rad'}; %... % relative flight-path angle ...
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%

alfa = '(alfamax*sin(alfatrig))';

% Elliptic cone
Cl = ['(Cl1*',alfa,'+ Cl0)'];
Cd = ['(Cd2*',alfa,'^2 + Cd1*',alfa,'+ Cd0)'];

% Quantities of Interest
rho = '(rho0*exp(-h/H))'; % Exponential Atmospheric Density [kg/m^3]
D   = ['(1/2*',rho,'*v^2*',Cd,'*Aref)']; % Drag Force [N]
L   = ['(1/2*',rho,'*v^2*',Cl,'*Aref)']; % Lift Force [N]

% Planar Equations of Motion
in.oc.stateRate = {'v*sin(gam)'; ...
				   'v*cos(gam)/(re+h)'; ...
				   ['-',D,'/mass - mu*sin(gam)/(re+h)^2']; ...
				   [L,'/(mass*v) + (v/(re+h) - mu/(v*(re+h)^2))*cos(gam)']};

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control = {'alfatrig','rad'}; % angle of attack control
in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximize/minimize
in.minimize = true;

% Path cost
in.oc.cost.path = {'1','s'};

% Terminal cost
in.oc.cost.terminal = {'0','s'};

% Initial cost
in.oc.cost.initial = {'0','s'};

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial constraint
in.oc.constraint.initial = {'h-x0(1)','m'; ...
	 						'thetta-x0(2)','rad'; ...
	 						'v-x0(3)','m/s'};

% Terminal constraint
in.oc.constraint.terminal = {'h-xf(1)','m'; ... 
                             'v-xf(3)','m/s'};

%%%%%%%%%%%%%%%
%% Constants %%
%%%%%%%%%%%%%%%

in.const.mu      = {42828.371901*1e9,'m^3/s^2'}; % Gravitational parameter, m^3/s^2
in.const.rho0    = {0.02,'kg/m^3'}; % Sea-level atmospheric density, kg/m^3
in.const.H       = {11100,'m'}; % Scale height for atmosphere of Earth, m
in.const.mass    = {92080,'kg'}; % Mass of vehicle, kg
in.const.re      = {3397000,'m'}; % Radius of planet, m
in.const.Aref    = {250,'m^2'}; % Reference area of vehicle, m^2
in.const.g0      = {3.71,'m/s^2'};
in.const.alfamax = {20*pi/180,'rad'};
in.const.Cl1     = {1.6756,'nd'};
in.const.Cl0     = {-0.2070,'nd'};
in.const.Cd2     = {2.04,'nd'};
in.const.Cd1     = {-0.3529,'nd'};
in.const.Cd0     = {0.0785,'nd'};
in.const.c1      = {10*pi/180,'rad'};
in.const.c0      = {0*pi/180,'rad'};
in.const.tol     = {1e-4,'nd'}; 
in.const.NMax    = {1e10,'nd'}; 
in.const.NMesh   = {500,'nd'};

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%

in.oc.guess.mode = 'auto';
in.oc.guess.timeIntegrate = 10;

% % Use automatic init
% Conditions at entry
in.oc.guess.initial.h      = 80e3;
in.oc.guess.initial.thetta = 0*pi/180;
in.oc.guess.initial.v      = 6000;
in.oc.guess.initial.gam    = -30*pi/180;

%%%%%%%%%%%%%%%%%%
%% Continuation %%
%%%%%%%%%%%%%%%%%%

in.cont.method = 1; % 1 = manually changing parameters

ind = 0;
%
ind = ind+1;
in.CONT{ind}.numCases = 100; % Number of steps in the continuation set
in.CONT{ind}.constraint.terminal.v = 4500;
in.CONT{ind}.constraint.terminal.h = 80e3;

%%%%%%%%%%%%%%%%%%%%%
% Continuation Set %% 
%%%%%%%%%%%%%%%%%%%%%

ind = ind+1;
in.CONT{ind}.numCases = 10;
in.CONT{ind}.const.tol = linspace(0,-(1e-4-1e-6),in.CONT{ind}.numCases);

return
