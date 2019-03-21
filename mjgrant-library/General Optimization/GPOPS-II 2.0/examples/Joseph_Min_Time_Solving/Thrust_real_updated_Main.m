%--------------- Thrust AoA Bang Bang ------------%
% Developed by: Kshitij Mall
% Last modified: 12/08/2015
% ------------------------------------------------%
clear; close all; clc

% Initial values
t0 = 0;
h0 = 20000;
thetta0 = 0*pi/180;
v0 = 1300;
gam0 = 0*pi/180;
mass0 = 1300;
aoa0 = 0*pi/180;
aoadot0 = 0*pi/180;
A0 = 0; %28;

% Final values
tf = 300;
hf = 0;
thettaf = 15.26*pi/180; %8.683/6378; %6.737*pi/180; %
vf = 3000;
gamf = 0*pi/180;
massf = 600;% 127005; % runs
aoaf = -10*pi/180;
aoadotf = 0;
Af = 0; %28;

% Minimum Values
tfMin = 0;
hMin = 0;
thettaMin = 0;
vMin = 1;
gamMin = -89.9*pi/180;
massMin = 600;
% aoaMin = -80*pi/180;  
aoaMin = -180*pi/180;  
aoadotMin = -5*pi/180;
AMin = 0; 

% Maximum Values
tfMax = 3000;
hMax = 380000;
thettaMax = 10*pi;
vMax = 10000;
gamMax = 89.9*pi/180;
massMax = mass0;
% aoaMax = 80*pi/180;
aoaMax = 180*pi/180;
aoadotMax = 5*pi/180;
AMax = 0.3;
%AMax = 1;

% Scramjet Vehicle:
% auxdata.cl   = [0.0049 0.496];               % Parameters for Lift Coefficient
% auxdata.cd   = [0.0007 0.0096 0.4747];        % Parameters for Drag Coefficient
% auxdata.cl   = [-0.0278 0.4639];               % Parameters for Lift Coefficient
% auxdata.cd   = [0.03 -0.0305 0.3216];        % Parameters for Drag Coefficient
auxdata.cl   = [0.1758 10.305];
auxdata.cd   = [0.26943 -0.4113 18.231];
auxdata.S    = 0.35;                    % Vehicle Reference Area (m^2) 
auxdata.rn = 1/12*0.3048; % Nose radius

% Earth
auxdata.Re   = 6378000;                     % Equatorial Radius of Earth (m)
auxdata.H    = 7500;                        % Density Scale Height (m)
auxdata.rho0 = 1.2;              % Sea Level Atmospheric Density (kg/m^3)
auxdata.mu   = 3.986e5*1e9;           % Earth Gravitational Parameter (m^3/s^2)
auxdata.g0 = 9.81; 

auxdata.Isp = 1600;
auxdata.Tmax = 1800000;

% New parameters

auxdata.y  = 1.4;
auxdata.R  = 287.058;
auxdata.T0  = 230;
auxdata.T4  = 1600;
auxdata.cp  = 1004;
auxdata.hpr  = 43903250;
auxdata.Mc  = 3;
auxdata.M1  = 4;
auxdata.M2  = 9.75;
auxdata.q1  = 40000;
auxdata.q2  = 150000;
auxdata.qscale = 1*10^3;
auxdata.T1  = 0;

%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [h0, thetta0, v0, gamMin, mass0, aoaMin];
bounds.phase.initialstate.upper = [h0, thetta0, v0, gamMax, mass0, aoaMax];
bounds.phase.state.lower = [hMin, thettaMin, vMin, gamMin, massMin,aoaMin];
bounds.phase.state.upper = [hMax, thettaMax, vMax, gamMax, massMax, aoaMax];
bounds.phase.finalstate.lower = [hf, thettaf, vMin, gamMin, massMin, aoaMin];
bounds.phase.finalstate.upper = [hf, thettaf, vMax, gamMax, massMax, aoaMax];
bounds.phase.control.lower = [aoadotMin, AMin];
bounds.phase.control.upper = [aoadotMax, AMax];
% bounds.phase.path.lower = 0.5; %50000;
% bounds.phase.path.upper = 900000; %100000;
guess.phase.time =  [t0; 90];
guess.phase.state(:,1) = [h0; h0];
guess.phase.state(:,2) = [thetta0; thetta0];
guess.phase.state(:,3) = [v0; v0];
guess.phase.state(:,4) = [0; gamMin];
guess.phase.state(:,5) = [mass0; massf];
guess.phase.state(:,6) = [aoa0; aoa0];
aoadotGuess         = [aoadotMin; aoadotMax];
AGuess           = [A0; AMin];
guess.phase.control = [aoadotGuess, AGuess];

    
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
setup.name                           = 'Thrust-real';
setup.functions.continuous           = @Thrust_real_updated_Continuous;
setup.functions.endpoint             = @Thrust_real_updated_Endpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 1;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'mumps';
%setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.supplier           = 'adigator';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';


%-------------------------------------------------------------------------%
%---------------------- Solve Problem using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
totaltime = tic;
output = gpops2(setup);
totaltime = toc(totaltime);