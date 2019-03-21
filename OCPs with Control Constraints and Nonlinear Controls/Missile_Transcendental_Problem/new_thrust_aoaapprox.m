function [in] = new_thrust_aoaapprox()
%
% This function creates input for the optimization problem,
% which is later used by other functions
%
% input : void
% output : in [structure]
% Developed by : J. Williams, K. Mall,  Dr. M. J. Grant and Thomas Antony
% Last modified: 29 March, 2017

%%%%%%%%%%%%%%%%%%%%%%%
%% Execution Control %%
%%%%%%%%%%%%%%%%%%%%%%%

in.oc.writeEquations = true; % Determine if we need to regenerate the equation files

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%
% These scaling parameters are altered automatically during continuation. This
% is just a first guess.
% Doubles must be equal to one (not dynamically updated during continuation)
in.autoScale = true;
in.bvpOptions = bvpset('AbsTol',1e-4,'RelTol',1e-4,'Nmax',1000000000000,'Stats','on');
% in.scale = {'m','x1'; ...
% 		  'rad',1; ...
% 		  's','x1/x3'; ...
% 		  'kg','x5'; ...
% 		  'nd',1}; % nd = nondimensional

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

%%%%%%%%%%%%
%% States %%
%%%%%%%%%%%%

% States
in.oc.state = {'h','m'; ... % altitude
		 'thetta','rad'; ... % longitude, positive eastward
		 'v','m/s'; ... % relative velocity
		 'gam','rad';... % relative flight-path angle 
         'mass','kg'};
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%

alfa = '(50*sin(alfatrig)*pi/180)';
% alfa = 'alfa';
% throttle = '((1+sin(Throttletrig))/2)';
% T = ['(Tmax*',throttle,')'];
% T = '(Tmax-Tmin)/2*sin(Throttletrig) + (Tmax+Tmin)/2';
% Elliptic cone
Cl = ['(0.4639*',alfa,'-0.0278)'];
Cd = ['(0.3216*',alfa,'^2 -0.0305*',alfa,'+0.03)'];
% % Elliptic cone
% Cl = ['(0.1*180/pi*',alfa,'-0)'];
% Cd = ['(0.0084*(180/pi)^2*',alfa,'^2 -0*',alfa,'+0.3)'];

% Quantities of Interest
rho = '(rho0*exp(-h/H))'; % Exponential Atmospheric Density [kg/m^3]
D = ['(1/2*',rho,'*v^2*',Cd,'*Aref)']; % Drag Force [N]
L = ['(1/2*',rho,'*v^2*',Cl,'*Aref)']; % Lift Force [N]

% Equations of Motion
%% Planar
in.oc.stateRate = {'v*sin(gam)'; ...
				'v*cos(gam)/(re+h)'; ...
                ['(Tmax-',D,')/mass - mu*sin(gam)/(re+h)^2']; ...
                ['(Tmax*',alfa, '+',L,')/(mass*v) + (v/(re+h) - mu/(v*(re+h)^2))*cos(gam)'];...
                '-Tmax/(g0*Isp)'};

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control = {'alfatrig','rad'};
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

% % Path cost
% in.oc.cost.path = {'0','m^2/s^2'};
% 
% % Terminal cost
% in.oc.cost.terminal = {'-v^2','m^2/s^2'};
% 
% % Initial cost
% in.oc.cost.initial = {'0','m^2/s^2'};
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial constraint
in.oc.constraint.initial = {'h-x0(1)','m'; ...
	 						'thetta-x0(2)','rad'; ...
	 						'v-x0(3)','m/s'; ...
                            'mass-x0(5)','kg'};

% Terminal constraint
in.oc.constraint.terminal = {'h-xf(1)','m'; ...
                             'thetta-xf(2)','rad'}; 

%%%%%%%%%%%%%%%%%%%%%%
%%     Constants    %%
%%%%%%%%%%%%%%%%%%%%%%

in.const.mu = {3.986e5*1e9,'m^3/s^2'}; % Gravitational parameter, m^3/s^2
in.const.rho0 = {1.2,'kg/m^3'}; % Sea-level atmospheric density, kg/m^3
in.const.H = {7500,'m'}; % Scale height for atmosphere of Earth, m
in.const.re = {6378000,'m'}; % Radius of planet, m
in.const.Aref = {557.4,'m^2'}; % Reference area of vehicle, m^2
in.const.Tmax  = {2000000,'kg*m/s^2'}; 
in.const.Isp  = {400,'s'}; % 1600
in.const.g0  = {9.80665,'m/s^2'};

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%
in.oc.guess.mode = 'auto';
in.oc.guess.timeIntegrate = 1;%10; 

% % Use automatic init
% Conditions at entry
in.oc.guess.initial.h = 26000; % 26000
in.oc.guess.initial.thetta = 0*pi/180;
in.oc.guess.initial.v = 1500; % 1600
in.oc.guess.initial.gam = 0*pi/180;
in.oc.guess.initial.mass = 136000;
in.oc.guess.costate = -0.1;%-0.1;

%%%%%%%%%%%%%%%%%%
%% Continuation %%
%%%%%%%%%%%%%%%%%%

in.cont.method = 1; % 1 = manually changing parameters

ind = 0;
%
ind = ind+1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set
in.CONT{ind}.constraint.terminal.h = 0;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 10; % Number of steps in the continuation set
in.CONT{ind}.constraint.terminal.thetta = 340/6378;  

return
