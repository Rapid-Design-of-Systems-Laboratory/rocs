%------------------------- Robot Arm Problem -----------------------------%
% This example is taken verbatim from the following reference:            %
%   Benchmarking Optimization Software with COPS Elizabeth D. Dolan       %
%   and Jorge J. More ARGONNE NATIONAL LABORATORY                         %
%-------------------------------------------------------------------------%
clear all
close all

%-------------------------------------------------------------------------%
%------------------ Provide Auxiliary Data for Problem -------------------%
%-------------------------------------------------------------------------%
auxdata.L = 5;

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0 = 0;
x10 = 4.5; x1f = 4.5;
x1min = 0; x1max = auxdata.L;
x20 = 0; x2f = 0;
x2min = -10*auxdata.L; x2max =  10*auxdata.L;
x30 = 0; x3f = 2*pi/3;
x3min = -pi; x3max =  pi;
x40 = 0; x4f = 0;
x4min = -50; x4max =  50;
x50 = pi/4; x5f = pi/4;
x5min =  0; x5max =  pi;
x60 = 0; x6f = 0;
x6min = -50; x6max =  50;
u1min = -1; u1max =  1;
u2min = -1; u2max =  1;
u3min = -1; u3max =  1;
tfmin = 0.1; tfmax= 10;

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfmin;
bounds.phase.finaltime.upper = tfmax;
bounds.phase.initialstate.lower = [x10, x20, x30, x40, x50, x60];
bounds.phase.initialstate.upper = [x10, x20, x30, x40, x50, x60];
bounds.phase.state.lower = [x1min, x2min, x3min, x4min, x5min, x6min];
bounds.phase.state.upper = [x1max, x2max, x3max, x4max, x5max, x6max];
bounds.phase.finalstate.lower = [x1f, x2f, x3f, x4f, x5f, x6f];
bounds.phase.finalstate.upper = [x1f, x2f, x3f, x4f, x5f, x6f];
bounds.phase.control.lower = [u1min, u2min, u3min];
bounds.phase.control.upper = [u1max, u2max, u3max];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess               = [t0; 5];
x1Guess               = [x10; x1f];
x2Guess               = [x20; x2f];
x3Guess               = [x30; x3f];
x4Guess               = [x40; x4f];
x5Guess               = [x50; x5f];
x6Guess               = [x60; x6f];
u1Guess               = [u1min; u1max];
u2Guess               = [u2min; u2max];
u3Guess               = [u3min; u3max];
guess.phase.state    = [x1Guess, x2Guess, x3Guess, x4Guess, x5Guess, x6Guess];
guess.phase.control  = [u1Guess, u2Guess, u3Guess];
guess.phase.time     = [tGuess];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-PattersonRao';
mesh.tolerance       = 1e-6;
mesh.maxiterations   = 10;
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 10;
N                    = 10;
N                    = 1;
mesh.phase.colpoints = 3*ones(1,N);
mesh.phase.fraction  = ones(1,N)/N;
mesh.phase.colpoints = 10*ones(1,N);
mesh.phase.fraction  = ones(1,N)/N;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.mesh                            = mesh;
setup.name                            = 'Robot-Arm-Problem';
setup.functions.endpoint              = @robotArmEndpoint;
setup.functions.continuous            = @robotArmContinuous;
setup.displaylevel                    = 2;
setup.auxdata                         = auxdata;
setup.bounds                          = bounds;
setup.guess                           = guess;
setup.nlp.solver                      = 'ipopt';
setup.nlp.ipoptoptions.linear_solver  = 'ma57';
setup.nlp.ipoptoptions.tolerance       = 1e-9;
setup.nlp.snoptoptions.tolerance       = 1e-8;
setup.derivatives.supplier            = 'sparseCD';
setup.derivatives.derivativelevel     = 'second';
setup.method                          = 'RPM-Differentiation';
setup.derivatives.derivativelevel     = 'second';
setup.mesh                            = mesh;
setup.scales.method                   = 'automatic-guess';

%-------------------------------------------------------------------------%
%----------------------- Solve Problem Using GPOPS2 ----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
