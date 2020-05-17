function [in] = main_heatrate_aoa()
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

in.oc.writeEquations = true; % Determine if we need to regenerate the equation files
in.useDeval          = true;

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%
% These scaling parameters are altered automatically during continuation. This
% is just a first guess.
% Doubles must be equal to one (not dynamically updated during continuation)
in.autoScale = true;
in.scale = {'m','x1'; ...
		    'rad',1; ...
		    's','x1/x3'; ...
            'kg','const.mass';...
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
in.oc.state = {'h','m'; ... % altitude
		       'thetta','rad'; ... % longitude, positive eastward
		       'v','m/s'; ... % relative velocity
		       'gam','rad'}; %... % relative flight-path angle ...
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%

rho = '(rho0*exp(-h/H))'; % Exponential Atmospheric Density [kg/m^3]
alfa = '(alfamax*sin(alfamctrig))';
Cl = ['(Cl1*',alfa,')'];
Cd = ['(Cd2*',alfa,'^2 + Cd0)'];

% Quantities of Interest
D = ['(1/2*',rho,'*v^2*',Cd,'*Aref)']; % Drag Force [N]
L = ['(1/2*',rho,'*v^2*',Cl,'*Aref)']; % Lift Force [N]

% Equations of Motion
%% Planar
in.oc.stateRate = {'v*sin(gam)'; ...
				   'v*cos(gam)/(re+h)'; ...
				   ['-',D,'/mass - mu*sin(gam)/(re+h)^2']; ...
                   [L,'/(mass*v) + (v/(re+h) - mu/(v*(re+h)^2))*cos(gam)']};

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control = {'alfamctrig','rad'}; % angle of attack trig control
in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path cost
in.oc.cost.path = {['epsilonh/cos(0.5*pi*k*sqrt(',rho,'/rn)*v^3/heatmax)'],'m/s^2'}; 

% Terminal cost 
in.oc.cost.terminal = {'-v','m/s'}; 

% Initial cost 
in.oc.cost.initial = {'0','m/s'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial constraint
in.oc.constraint.initial = {'h-x0(1)','m'; ...
	 						'thetta-x0(2)','rad'; ...
	 						'v-x0(3)','m/s';...
                            'gam-x0(4)','rad'; 
 						    };

% Terminal constraint
in.oc.constraint.terminal = {'h-xf(1)','m'; ... 
                             'thetta-xf(2)','rad';
							};

%%%%%%%%%%%%%%%%%%%%%%
%%     Constants    %%
%%%%%%%%%%%%%%%%%%%%%%
in.const.mu       = {3.986e5*1e9,'m^3/s^2'}; % Gravitational parameter
in.const.rho0     = {1.2,'kg/m^3'}; % Sea-level atmospheric density
in.const.H        = {7500,'m'}; % Scale height for atmosphere of Earth
in.const.mass     = {750/2.2046226,'kg'}; % Mass of vehicle
in.const.re       = {6.3781e6,'m'}; % Radius of planet
in.const.Aref     = {pi*(24*.0254/2)^2,'m^2'}; % Reference area of vehicle
in.const.g        = {9.80665,'m/s^2'};
in.const.gmax     = {150,'nd'}; 
in.const.Cl1      = {1.5658,'nd'};
in.const.Cd0      = {0.0612,'nd'};
in.const.Cd2      = {1.6537,'nd'};
in.const.alfamax  = {80*pi/180,'rad'};
in.const.heatmax  = {10000e4,'kg/s^3'}; 
in.const.k        = {1.74153e-4,'sqrt(kg)/m'}; % heat rate coefficient
in.const.rn       = {1/12*0.3048,'m'}; % Nose radius
in.const.epsilonh = {1,'m/s^2'}; 
in.const.tol      = {1e-6,'nd'}; 
in.const.NMax     = {1e13,'nd'}; 
in.const.NMesh    = {500,'nd'};

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%
in.oc.guess.mode = 'auto';
in.oc.guess.timeIntegrate = 1;

% % Use automatic init
% Conditions at entry
in.oc.guess.initial.h = 80e3;
in.oc.guess.initial.thetta = 0*pi/180;
in.oc.guess.initial.v = 5000;
in.oc.guess.initial.gam = -89*pi/180 ;
in.oc.costates = -0.1;

%%%%%%%%%%%%%%%%%%
%% Continuation %%
%%%%%%%%%%%%%%%%%%

in.cont.method = 1; % 1 = manually changing parameters

ind = 0;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 1
%%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.CONT{ind}.numCases = 10; 
in.CONT{ind}.constraint.terminal.h = 15e3;
in.CONT{ind}.constraint.terminal.thetta = 0.01*pi/180;
% in.CONT{ind}.constraint.initial.h = 50e3;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 2
%%%%%%%%%%%%%%%%%%%%%%
ind = ind + 1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set
in.CONT{ind}.constraint.initial.gam = -60*pi/180;
in.CONT{ind}.constraint.terminal.thetta = 0.5*pi/180;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 3
%%%%%%%%%%%%%%%%%%%%%%
ind = ind + 1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 10;
in.CONT{ind}.constraint.terminal.thetta = 1*pi/180;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 4
%%%%%%%%%%%%%%%%%%%%%%

ind = ind + 1;
in.CONT{ind}.numCases = 50; % Number of steps in the continuation set
in.cont.method(ind) = 1;
in.CONT{ind}.const.heatmax = linspace(0,-8000e4,in.CONT{ind}.numCases);
% in.CONT{ind}.const.epsilonh = linspace(0,-0.09999,in.CONT{ind}.numCases);

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 5
%%%%%%%%%%%%%%%%%%%%%%
ind = ind + 1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set % 44
in.cont.method(ind) = 1;
in.CONT{ind}.const.alfamax = linspace(0,-40*pi/180,in.CONT{ind}.numCases);

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 6
%%%%%%%%%%%%%%%%%%%%%%

ind = ind + 1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set
in.cont.method(ind) = 1;
in.CONT{ind}.const.epsilonh = linspace(0,-(1-0.01),in.CONT{ind}.numCases);

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 7
%%%%%%%%%%%%%%%%%%%%%%

ind = ind + 1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set
in.cont.method(ind) = 1;
in.CONT{ind}.const.epsilonh = linspace(0,-(0.01-0.001),in.CONT{ind}.numCases);

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 8
%%%%%%%%%%%%%%%%%%%%%%

ind = ind + 1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set
in.cont.method(ind) = 1;
in.CONT{ind}.const.epsilonh = linspace(0,-(0.001-0.0001),in.CONT{ind}.numCases);

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %% 9
%%%%%%%%%%%%%%%%%%%%%%

ind = ind + 1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set
in.cont.method(ind) = 1;
in.CONT{ind}.const.epsilonh = linspace(0,-(0.0001-0.00001),in.CONT{ind}.numCases);
return
