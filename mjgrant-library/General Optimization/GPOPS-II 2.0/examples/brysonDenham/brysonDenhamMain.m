%------------------------- Bryson-Denham Problem -------------------------%
% This problem is taken from the following reference:                     %
% Bryson, A. E. and Ho, Y-C, "Applied Optimal Control:  Optimization,     %
% Estimation, and Control," Hemisphere Publishing, 1975.                  %
%-------------------------------------------------------------------------%
clear all
close all
clc

%-------------------------------------------------------------------------%
%------------- Provide and Set Up All Bounds for Problem -----------------%
%-------------------------------------------------------------------------%
L                               = 1/9;
t0                              = 0;     
tf                              = 1;
x0                              = 0;
xf                              = 0;
xmin                            = 0;
xmax                            = L;
v0                              = 1;
vf                              = -v0;
vmin                            = -10;
vmax                            = 10;
umin                            = -10;
umax                            = 5;
bounds.phase.initialtime.lower  = t0;
bounds.phase.initialtime.upper  = t0;
bounds.phase.finaltime.lower    = tf;
bounds.phase.finaltime.upper    = tf;
bounds.phase.initialstate.lower = [x0,v0];
bounds.phase.initialstate.upper = [x0,v0];
% bounds.phase.state.lower        = [xmin,vmin];
% bounds.phase.state.upper        = [xmax,vmax];
bounds.phase.state.lower        = [xmin,vmin];
bounds.phase.state.upper        = [10,vmax];
bounds.phase.finalstate.lower   = [xf,vf];
bounds.phase.finalstate.upper   = [xf,vf];
bounds.phase.control.lower      = [umin];
bounds.phase.control.upper      = [umax];
bounds.phase.integral.lower     = 0;
bounds.phase.integral.upper     = 10;
bounds.phase.path.lower         = [xmin];
bounds.phase.path.upper         = [xmax];

t1 = (0:0.3/100:0.3).';
t2 = (0.3:0.4/100:0.7).';
t3 = (0.7:0.3/100:1).';

t1 = t1(1:end-1);
t2 = t2(1:end-1);
t3 = t3(1:end-1);

xE{1}(:,1) = L*(1-(1-t1/(3*L)).^3);
xE{2}(:,1) = L*ones(size(t2));
xE{3}(:,1) = L*(1-(1-(1-t3)/(3*L)).^3);
xE{1}(:,2) = (1-t1/(3*L)).^2;
xE{2}(:,2) = zeros(size(t2));
xE{3}(:,2) = -(1-(1-t3)/(3*L)).^2;
uE{1}(:,1) = -(2/(3*L))*(1-t1/(3*L));
uE{2}(:,1) = zeros(size(t2));
uE{3}(:,1) = -(2/(3*L))*(1-(1-t3)/(3*L));

tg = [t1; t2; t3];
xg = [xE{1}; xE{2}; xE{3}];
ug = [uE{1}; uE{2}; uE{3}];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess               = [t0;tf];
xGuess               = [x0;xf];
vGuess               = [v0;vf];
uGuess               = [umin;0];
guess.phase.time     = tGuess;
guess.phase.state    = [xGuess,vGuess];
guess.phase.control  = uGuess;
guess.phase.integral = 0;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method            = 'hp-PattersonRao';
mesh.tolerance         = 1e-10;
mesh.maxiterations     = 20;
mesh.colpointsmin      = 3;
mesh.colpointsmax      = 16;
NumIntervals           = 4;
mesh.phase.colpoints   = 6*ones(1,NumIntervals);
mesh.phase.fraction    = ones(1,NumIntervals)/NumIntervals;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                           = 'Bryson-Denham-Problem';
setup.functions.continuous           = @brysonDenhamContinuous;
setup.functions.endpoint             = @brysonDenhamEndpoint;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.tolerance     = 1e-8;
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.method                         = 'RPM-Differentiation';
setup.mesh                           = mesh;

%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
