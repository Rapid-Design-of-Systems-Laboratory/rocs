%--------------- Aircraft Noise ------------%
% Developed by: Kshitij Mall and Janav Udani
% Last modified: 09/04/2016
% ------------------------------------------%
clear; close all; clc

% Initial values
t0 = 0;
x0 = 0; 
z0 = 1197; 
v0 = 124; 
gam0 = 0*pi/180;
aoa0 = 0*pi/180;
T0 = 3420;

% Final values
tf = 180;
xf = 15000; 
zf = 0;
vf = 77.5; 
gamf = 0*pi/180;
aoaf = 0*pi/180;
Tf = 300; 

% Minimum Values
tfMin = 0;
xMin = 0;
zMin = 0;
vMin = 71.5;
gamMin = -90*pi/180;
aoaMin = -10*pi/180;
TMin = 300; 

% Maximum Values
tfMax = 1000;
xMax = 80000;
zMax = 10000;
vMax = 282;
gamMax = 90*pi/180; 
aoaMax = 10*pi/180; 
TMax = 3420;

% Aircraft
auxdata.cl   = [0.1667 6.231];%  Parameters for Lift Coefficient
auxdata.cd   = [0.07351 -0.08617 1.996]; % Parameters for Drag Coefficient
auxdata.S    = 18.2; % Vehicle Reference Area (m^2)
auxdata.mass = 7180/9.81; % Mass

% Earth
auxdata.H    = 7500; % Density Scale Height (m)
auxdata.rho0 = 1.2; % Sea Level Atmospheric Density (kg/m^3)
auxdata.g0 = 9.81;
auxdata.C3 = 18.73;
auxdata.C1 = 0.226;
auxdata.C2 = 5.2e6;

%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [x0, z0, v0, gam0];
bounds.phase.initialstate.upper = [x0, z0, v0, gam0];
bounds.phase.state.lower = [xMin, zMin, vMin, gamMin];
bounds.phase.state.upper = [xMax, zMax, vMax, gamMax];
bounds.phase.finalstate.lower = [xf, zf, vf, gamf];
bounds.phase.finalstate.upper = [xf, zf, vf, gamf];
bounds.phase.control.lower = [aoaMin, TMin];
bounds.phase.control.upper = [aoaMax, TMax];
bounds.phase.integral.lower     = 10^1;
bounds.phase.integral.upper     = 10^25;
guess.phase.time =  [t0; tf];
guess.phase.state(:,1) = [x0; xf];

guess.phase.state(:,2) = [z0; zf];
guess.phase.state(:,3) = [v0; vf];
guess.phase.state(:,4) = [gam0; gamf];
aoaGuess         = [aoaMin; aoaMax];
TGuess           = [TMax; TMin];
guess.phase.control = [aoaGuess, TGuess];
guess.phase.integral = 10^12;
    
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-PattersonRao';
mesh.maxiterations = 100;
mesh.colpointsmin = 3;
mesh.colpointsmax = 10;
mesh.tolerance    = 1e-6;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Aircraft_Noise_2DOF';
setup.functions.continuous           = @Aircraft_Noise_Continuous_2DOF_Pop;
setup.functions.endpoint             = @Aircraft_Noise_Endpoint_2DOF_Pop;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 1;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'mumps';
% setup.nlp.solver                     = 'snopt';
% setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
setup.derivatives.supplier           = 'adigator';
setup.derivatives.method             = 'adigator';

%-------------------------------------------------------------------------%
%---------------------- Solve Problem using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
totaltime = tic;
output = gpops2(setup);
totaltime = toc(totaltime);