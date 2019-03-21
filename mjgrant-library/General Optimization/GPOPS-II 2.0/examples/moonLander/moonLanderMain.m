%----------------------- Moon-Lander Problem -----------------------------%
% This example can be found in the following reference:                   %
% Meditch, J., "On the Problem of Optimal Thrust Programming for a Soft   %
% Lunar Landing," IEEE Transactions on Automatic Control, Vol. 9,% No. 4, %
% 1964, pp. 477-484.                                                      %
%-------------------------------------------------------------------------%
clear all; close all; clc

%-------------------------------------------------------------------------%
%------------------ Provide Auxiliary Data for Problem -------------------%
%-------------------------------------------------------------------------%
auxdata.g = 1.5;

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0min = 0;  t0max = 0;
tfmin = 0;  tfmax = 200;
h0 = 10;    hf = 0;
v0 = -2;    vf = 0;
hmin = 0;   hmax =  20;
vmin = -10; vmax =  10;
umin = 0;   umax = 3;

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0min;
bounds.phase.initialtime.upper = t0max;
bounds.phase.finaltime.lower = tfmin;
bounds.phase.finaltime.upper = tfmax;
bounds.phase.initialstate.lower = [h0, v0];
bounds.phase.initialstate.upper = [h0, v0];
bounds.phase.state.lower = [hmin, vmin];
bounds.phase.state.upper = [hmax, vmax];
bounds.phase.finalstate.lower = [hf, vf];
bounds.phase.finalstate.upper = [hf, vf];
bounds.phase.control.lower = [umin];
bounds.phase.control.upper = [umax];
bounds.phase.integral.lower = [-100];
bounds.phase.integral.upper = [100];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess               = [t0min; 5];
hGuess               = [h0; hf];
vGuess               = [v0; vf];
uGuess               = [umin; umin];
guess.phase.state    = [hGuess, vGuess];
guess.phase.control  = [uGuess];
guess.phase.time     = [tGuess];
guess.phase.integral = 10;

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
setup.name                        = 'Soft-Lunar-Landing';
setup.functions.endpoint          = @moonLanderEndpoint;
setup.functions.continuous        = @moonLanderContinuous;
setup.displaylevel                = 2;
setup.auxdata                     = auxdata;
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
