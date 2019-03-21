%------------------------- Bryson-Denham Problem -------------------------%
% This problem is taken from the following reference:                     %
% Bryson, A. E. and Ho, Y-C, "Applied Optimal Control:  Optimization,     %
% Estimation, and Control," Hemisphere Publishing, 1975.                  %
%-------------------------------------------------------------------------%
clear all; clc

%--------------------------------------------------------------------------%
%---------------- Set Up Auxiliary Data for Problem -----------------------%
%--------------------------------------------------------------------------%
auxdata.T  = 0.1405;
auxdata.m0 = 1;
auxdata.dm = 0.0749;
auxdata.mu = 1;

%--------------------------------------------------------------------------%
%--------------- Set Up Bounds on State, Control, and Time ----------------%
%--------------------------------------------------------------------------%
t0          = 0; 
tf          = 3.32;
r0          = 1; 
theta0      = 0; 
vr0         = 0; 
vtheta0     = 1; 
vrf         = 0;
rmin        = 1; 
rmax        = 10;
thetamin    = -pi; 
thetamax    = pi;
vrmin       = -10; 
vrmax       = 10;
vthetamin   = -pi; 
vthetamax   = pi;
u1min       = -1.3; 
u1max       = 1.3;
u2min       = -1.3; 
u2max       = 1.3;
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = [r0, theta0, vr0, vtheta0];
bounds.phase.initialstate.upper = [r0, theta0, vr0, vtheta0];
bounds.phase.state.lower = [rmin, thetamin, vrmin, vthetamin];
bounds.phase.state.upper = [rmax, thetamax, vrmax, vthetamax];
bounds.phase.finalstate.lower = [rmin, thetamin, vrf, vthetamin];
bounds.phase.finalstate.upper = [rmax, thetamax, vrf, vthetamax];
bounds.phase.control.lower = [u1min, u2min];
bounds.phase.control.upper = [u1max, u2max];
bounds.phase.path.lower  = 1;
bounds.phase.path.upper  = 1;
bounds.eventgroup.lower = [0];
bounds.eventgroup.upper = [0];

%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Guess ---------------------------%
%--------------------------------------------------------------------------%
tGuess      = [t0; tf];
rGuess      = [r0; 1.5];
thetaGuess  = [theta0; pi];
vrGuess     = [vr0; vrf];
vthetaGuess = [vtheta0; .5];
u1Guess   = [0; -1];
u2Guess   = [1; 0];
guess.phase.time    = [tGuess];
guess.phase.state   = [rGuess, thetaGuess, vrGuess, vthetaGuess];
guess.phase.control = [u1Guess, u2Guess];

%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Mesh ----------------------------%
%--------------------------------------------------------------------------%
N = 10;
meshphase.colpoints = 4*ones(1,N);
meshphase.fraction   = ones(1,N)/N;

%--------------------------------------------------------------------------%
%-------------------------- Set Up for Solver -----------------------------%
%--------------------------------------------------------------------------%
setup.name = 'Orbit-Raising-Problem';
setup.functions.continuous = @orbitRaisingContinuous;
setup.functions.endpoint = @orbitRaisingEndpoint;

setup.displaylevel = 2;
setup.nlp.solver = 'ipopt';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;

setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';

setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-6;
setup.mesh.phase = meshphase;
setup.scales.method = 'automatic-hybrid';

%--------------------------------------------------------------------------%
%-------------------- Solve Problem and Extract Solution ------------------%
%--------------------------------------------------------------------------%
output = gpops2(setup);
