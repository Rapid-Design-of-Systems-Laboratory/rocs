%---------------------- Double-Integrator Problem ------------------------%
% This example can be found in any good book on optimal control theory    %
%-------------------------------------------------------------------------%
clear all; close all; clc

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0min = 0;  t0max = 0;
tfmin = 0;  tfmax = 200;
x10   = +2; x1f   = 0;
x20   = +1; x2f   = 0;
x1min = -20; x1max =  20;
x2min = -20; x2max =  20;
umin = -1;   umax = +1;

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0min;
bounds.phase.initialtime.upper = t0max;
bounds.phase.finaltime.lower = tfmin;
bounds.phase.finaltime.upper = tfmax;
bounds.phase.initialstate.lower = [x10, x20];
bounds.phase.initialstate.upper = [x10, x20];
bounds.phase.state.lower = [x1min, x2min];
bounds.phase.state.upper = [x1max, x2max];
bounds.phase.finalstate.lower = [x1f, x1f];
bounds.phase.finalstate.upper = [x1f, x1f];
bounds.phase.control.lower = [umin];
bounds.phase.control.upper = [umax];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess               = [0; 5];
x1Guess              = [x10; x1f];
x2Guess              = [x20; x2f];
uGuess               = [umin; umin];
guess.phase.state    = [x1Guess, x2Guess];
guess.phase.control  = [uGuess];
guess.phase.time     = [tGuess];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-PattersonRao';
mesh.tolerance       = 1e-6;
mesh.maxiterations    = 20;
mesh.colpointsmin    = 4;
mesh.colpointsmax    = 10;
mesh.phase.colpoints = 4;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.mesh                        = mesh;
setup.name                        = 'Double-Integrator-Min-Time';
setup.functions.endpoint          = @doubleIntegratorEndpoint;
setup.functions.continuous        = @doubleIntegratorContinuous;
setup.displaylevel                = 2;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.method                      = 'RPM-Differentiation';

%-------------------------------------------------------------------------%
%----------------------- Solve Problem Using GPOPS2 ----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
