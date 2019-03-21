%-------- Minimum-Time Reorientation of an Asymmetric Rigid Body ---------%
% This problem is taken from the following reference:                     %
% Betts, John T. "Practical Methods for Optimal Control and Estimation    %
% Using Nonlinear Programming", 2nd Edition, Page 299-304                 %
%-------------------------------------------------------------------------%
clear all; clc

%--------------------------------------------------------------------------%
%---------------- Set Up Auxiliary Data for Problem -----------------------%
%--------------------------------------------------------------------------%
auxdata.Ix = 5621; % Kg*m^2
auxdata.Iy = 4547; % Kg*m^2
auxdata.Iz = 2364; % Kg*m^2
auxdata.ang = 150*pi/180; % 150 degree maneuver about X-axis

%--------------------------------------------------------------------------%
%--------------- Set Up Bounds on State, Control, and Time ----------------%
%--------------------------------------------------------------------------%
% Time
t0          = 0; 
tfMin       = 0;
tfMax       = 30;
% Initial States
q10         = 0; 
q20         = 0; 
q30         = 0; 
w10         = 0; 
w20         = 0; 
w30         = 0;
% Final States
q1f         = sin(auxdata.ang/2); 
q2f         = 0; 
q3f         = 0; 
w1f         = 0; 
w2f         = 0; 
w3f         = 0; 
% Initial & Final Control
q40         = 1;
q4f         = cos(auxdata.ang/2); 
% Minimum & Max States
q1min       = -1;       q1max        = 1;
q2min       = -1;       q2max        = 1;
q3min       = -1;       q3max        = 1;
w1min       = -1;       w1max        = 1;
w2min       = -1;       w2max        = 1;
w3min       = -1;       w3max        = 1;
% Minimum & Max Controls
q4min       = -1;       q4max        = 1;
u1min       = -50;      u1max        = 50;
u2min       = -50;      u2max        = 50;
u3min       = -50;      u3max        = 50;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [q10, q20, q30, w10, w20, w30];
bounds.phase.initialstate.upper = [q10, q20, q30, w10, w20, w30];
bounds.phase.state.lower = [q1min, q2min, q3min, w1min, w2min, w3min]; 
bounds.phase.state.upper = [q1max, q2max, q3max, w1max, w2max, w3max]; 
bounds.phase.finalstate.lower = [q1f, q2f, q3f, w1f, w2f, w3f];
bounds.phase.finalstate.upper = [q1f, q2f, q3f, w1f, w2f, w3f];
bounds.phase.control.lower = [q4min, u1min, u2min, u3min];
bounds.phase.control.upper = [q4max, u1max, u2max, u3max];
bounds.phase.path.lower  = 1;
bounds.phase.path.upper  = 1;

%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Guess ---------------------------%
%--------------------------------------------------------------------------%
% Time
tGuess      = [t0; tfMax];
% States
q1Guess     = [q10; q1f];
q2Guess     = [q20; q2f];
q3Guess     = [q30; q3f];
w1Guess     = [w10; w1f];
w2Guess     = [w20; w2f];
w3Guess     = [w30; w3f];
% Controls
q4Guess     = [q40; q4f];
u1Guess     = [ 50; 0];
u2Guess     = [-50; 0];
u3Guess     = [ 50; 0];

guess.phase.time    = [tGuess];
guess.phase.state   = [q1Guess, q2Guess, q3Guess, w1Guess, w2Guess, w3Guess];
guess.phase.control = [q4Guess, u1Guess, u2Guess, u3Guess];

%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Mesh ----------------------------%
%--------------------------------------------------------------------------%
N = 10;
mesh.method = 'hp-PattersonRao';
mesh.phase.colpoints = 4*ones(1,N);
mesh.phase.fraction   = ones(1,N)/N;
mesh.tolerance = 1e-8;

%--------------------------------------------------------------------------%
%-------------------------- Set Up for Solver -----------------------------%
%--------------------------------------------------------------------------%
setup.name = 'Reorientation-Problem';
setup.functions.continuous = @reorientationContinuous;
setup.functions.endpoint = @reorientationEndpoint;
setup.method = 'RPM-Differentiation';
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.mesh = mesh;

%--------------------------------------------------------------------------%
%-------------------- Solve Problem and Extract Solution ------------------%
%--------------------------------------------------------------------------%
output = gpops2(setup);

