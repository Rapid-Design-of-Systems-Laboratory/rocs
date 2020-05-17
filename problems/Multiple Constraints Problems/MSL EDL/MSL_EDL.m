function [in] = MSL_EDL()
%
% This function creates input for the optimization problem,
% which is later used by other functions
% input : void
% output : in [structure]
% Developed by : Dr. K. Mall 
% Last modified: Jan 28, 2020

%%%%%%%%%%%%%%%%%%%%%%%
%% Execution Control %%
%%%%%%%%%%%%%%%%%%%%%%%

in.oc.writeEquations = false; % Determine if we need to regenerate the equation files

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%

in.autoScale = true;
% in.scale = {'m','x1'; ...
%             'rad',1; ...
%             's','x1/x2'; ...
%             'kg','const.mass'; ...
%             'nd',1}; % nd = nondimensional

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

in.oc.state = {'h','m'; ... % Radial position magnitude
               'v','m/s'; ... % relative velocity
               'gam','rad'}; % Relative flight-path angle

%%%%%%%%%%%%%%%%%%%%%%%%
%% Vehicle Parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%
rho = '(rho0*exp(-h/H))'; % Exponential atmospheric density

% Find vehicle aerodynamic parameters
D = ['(1/2*',rho,'*v^2*Cd*Aref)']; % Drag Force
L = ['(1/2*',rho,'*v^2*Cd*LbyD*Aref)']; % Lift Force

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D
Ft = D; % Force along velocity vector
Fn = L; % Force perpendicular to velocity vector
control = '(0.683*sin(MSLctrl) + 0.183)';

in.oc.stateRate = {'v*sin(gam)'; ...
                   ['-',Ft,'/mass - muu*sin(gam)/(rm + h)^2']; ...
                   [Fn,'*',control,'/(mass*v) - muu/(v*(rm+h)^2)*cos(gam) + v/(rm+h)*cos(gam)']};

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control             = {'MSLctrl','rad'}; % angle of attack control
in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximize terminal altitude %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path cost
gload = ['(sqrt(',L,'^2 + ',D,'^2)/(mass*g))'];
q     = ['(0.5*',rho,'*v^2)'];
Q     = ['(k*sqrt(',rho,'/rn)*v^3)'];

% % For control constraint only
% in.oc.cost.path     = {'epsilon*cos(MSLctrl)','m/s'};

% For control and state path constraint
in.oc.cost.path     = {['epsilon*cos(MSLctrl) + epsq*sec(0.5*pi*(',q,')/qmax) + epsg*sec(0.5*pi*(',gload,')/gmax) + epsQ*sec(0.5*pi*(',Q,')/Qmax)'],'m/s'};

% Terminal cost
in.oc.cost.terminal = {'-h','m'};

% Initial cost
in.oc.cost.initial  = {'0','m'};

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial constraint
in.oc.constraint.initial = {'h-x0(1)','m'; ...
                            'v-x0(2)','m/s';
                            'gam-x0(3)','rad'};

% Terminal constraint
in.oc.constraint.terminal = {'v-xf(2)','m/s'};
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      HMME Mission Parameters    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants

in.const.rm      = {3397000,'m'}; % b
in.const.muu     = {4.284*1e13,'m^3/s^2'}; % Gravitational parameter
in.const.rho0    = {0.0158,'kg/m^3'}; % Sea-level atmospheric density
in.const.H       = {9354,'m'}; % Scale height for atmosphere of Earth
in.const.mass    = {3300,'kg'}; % Mass of vehicle
in.const.Aref    = {15.9,'m^2'}; % Reference area of vehicle
in.const.rn      = {0.6,'m'}; % Nose radius 
in.const.k       = {1.9027e-4,'sqrt(kg)/m^2'}; % heat rate coefficient
in.const.epsilon = {1,'m/s'}; % scaling factor for smoothed bang-bang control
in.const.epsQ    = {1,'m/s'}; % scaling factor for smoothed bang-bang control
in.const.epsg    = {1,'m/s'}; % scaling factor for smoothed bang-bang control
in.const.epsq    = {1,'m/s'}; % scaling factor for smoothed bang-bang control
in.const.LbyD    = {0.24, 'nd'};
in.const.Cd      = {1.45, 'nd'};
in.const.gmax    = {50, 'nd'};
in.const.Qmax    = {200*1e4, 'kg/(s^3*m)'}; 
in.const.qmax    = {100e3, 'kg/(m*s^2)'};
in.const.g       = {9.81, 'm/s^2'};
in.const.umin    = {cosd(120), 'nd'};
in.const.umax    = {cosd(30), 'nd'};
in.const.tol     = {1e-4,'nd'}; 
in.const.NMax    = {1e13,'nd'}; 

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%

% in.oc.initialGuessFunc = @getInitialGuess;
in.oc.guess.mode          = 'auto';
in.oc.guess.timeIntegrate = 10; 

% Use automatic init
% Conditions at entry
in.oc.guess.initial.h      = 50000;
in.oc.guess.initial.v      = 6000;
in.oc.guess.initial.gam    = -11.5*pi/180;
in.oc.costates             = -0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Execution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in.cont.method = 1; 
ind = 0;
    
%%%%%%%%%%%%%%%%%%%%%
% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.CONT{ind}.numCases = 50;
in.CONT{ind}.constraint.terminal.v = 540; 
in.CONT{ind}.constraint.initial.h = 125000; 

%%%%%%%%%%%%%%%%%%%%%
% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.CONT{ind}.numCases = 100;
in.CONT{ind}.const.Qmax = linspace(0,-(200e4-70e4),in.CONT{ind}.numCases);
% 
%%%%%%%%%%%%%%%%%%%%%
% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.CONT{ind}.numCases = 100;
in.CONT{ind}.const.gmax = linspace(0,-(50-5),in.CONT{ind}.numCases);
% 
%%%%%%%%%%%%%%%%%%%%%
% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.CONT{ind}.numCases = 100;
in.CONT{ind}.const.qmax = linspace(0,-(100e3-10e3),in.CONT{ind}.numCases);

%%%%%%%%%%%%%%%%%%%%%
% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.CONT{ind}.numCases = 300;
in.CONT{ind}.const.epsilon = linspace(0,-(1-1e-4),in.CONT{ind}.numCases);
in.CONT{ind}.const.epsQ = linspace(0,-(1-1e-2),in.CONT{ind}.numCases);
in.CONT{ind}.const.epsg = linspace(0,-(1-1e-2),in.CONT{ind}.numCases);
in.CONT{ind}.const.epsq = linspace(0,-(1-1e-2),in.CONT{ind}.numCases);

%%%%%%%%%%%%%%%%%%%%%
% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%
ind = ind+1;
in.CONT{ind}.numCases = 1000;
in.CONT{ind}.const.epsilon = linspace(0,-(1e-4-1e-6),in.CONT{ind}.numCases);
in.CONT{ind}.const.epsQ = linspace(0,-(1e-2-1e-6),in.CONT{ind}.numCases);
in.CONT{ind}.const.epsg = linspace(0,-(1e-2-1e-6),in.CONT{ind}.numCases);
in.CONT{ind}.const.epsq = linspace(0,-(1e-2-1e-6),in.CONT{ind}.numCases);
in.CONT{ind}.const.tol = linspace(0,-(1e-4-1e-6),in.CONT{ind}.numCases);

return
