%----------------------- Hyper-Sensitive Problem -------------------------%
% This example is taken from the following reference:                     %
% Rao, A. V., and Mease, K. D., "Eigenvector Approximate Dichotomic Basis %
% Methods for Solving Hyper-Sensitive Optimal Control Problems," Optimal  %
% Control Applications and Methods, Vol. 21, No. 1., January-February,    %
% 2000, pp. 1-17.                                                         %
%-------------------------------------------------------------------------%
clear all; clc

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0 = 0;
tf = 10000;
x0 = 1.5;
xf = 1; 
xMin = -50;
xMax = +50;
uMin = -50;
uMax = +50;

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = x0; 
bounds.phase.initialstate.upper = x0;
bounds.phase.state.lower = xMin;
bounds.phase.state.upper = xMax;
bounds.phase.finalstate.lower = xf; 
bounds.phase.finalstate.upper = xf;
bounds.phase.control.lower = uMin; 
bounds.phase.control.upper = uMax;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 100000;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time     = [t0; tf]; 
guess.phase.state    = [x0; xf];
guess.phase.control  = [0; 0];
guess.phase.integral = 0;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-7;
mesh.maxiterations   = 45;
mesh.colpointsmin    = 2;
mesh.colpointsmax    = 14;
mesh.phase.colpoints = 4*ones(1,10);
mesh.phase.fraction  = 0.1*ones(1,10);

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                           = 'Hyper-Sensitive-Problem';
setup.functions.continuous           = @hyperSensitiveContinuous;
setup.functions.endpoint             = @hyperSensitiveEndpoint;
setup.displaylevel                   = 2;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.nlp.solver                     = 'snopt';
setup.nlp.snoptoptions.tolerance     = 1e-10;
setup.nlp.snoptoptions.maxiterations = 20000;
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-10;
setup.derivatives.supplier           = 'adigator';
setup.derivatives.derivativelevel    = 'second';
setup.method                         = 'RPM-Differentiation';

%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
toc
