%---------------- Tumor Anti-Angiogenesis Problem --------------------%
% This example is taken from the following reference:                 %
% Ledzewicz, U. and Schattler, H, "Analysis of Optimal Controls for a %
% Mathematical Model of Tumour Anti-angiogenesis," Optimal Control    %
% Applications and Methods, Vol. 29, 2008, pp. 41-57.                 %
%---------------------------------------------------------------------%
clear all
clc

% Parameters:
%-------------------------------------------------------------------%
%-------------------- Data Required by Problem ---------------------%
%-------------------------------------------------------------------%
auxdata.zeta = 0.084;       % per day
auxdata.b = 5.85;           % per day
auxdata.d = 0.00873;        % per mm^2 per day
auxdata.G = 0.15;           % per mg of dose per day
auxdata.mu = 0.02;          % per day
a = 75;           
A = 15;

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
pMax = ((auxdata.b-auxdata.mu)/auxdata.d)^(3/2);
pMin = 0.1;
qMax = pMax;
qMin = pMin;
yMax = A;
yMin = 0;
uMax = a;
uMin = 0;
t0Max = 0;
t0Min = 0;
tfMax = 5;
tfMin = 0.1;
p0 = pMax/2;
q0 = qMax/4;
y0 = 0;

bounds.phase.initialtime.lower = t0Min;
bounds.phase.initialtime.upper = t0Max;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [p0, q0];
bounds.phase.initialstate.upper = [p0, q0];
bounds.phase.state.lower = [pMin, qMin];
bounds.phase.state.upper = [pMax, qMax];
bounds.phase.finalstate.lower = [pMin, qMin];
bounds.phase.finalstate.upper = [pMax, qMax];
bounds.phase.control.lower = uMin;
bounds.phase.control.upper = uMax;
bounds.phase.integral.lower = [0];
bounds.phase.integral.upper = [A];
guess.phase.time    = [0; 1];
guess.phase.state   = [[p0; pMax],[q0; qMax]];
guess.phase.control = [uMax; uMax];
guess.phase.integral = [7.5];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-6;
mesh.phase.colpoints = 4*ones(1,10);
mesh.phase.fraction  = 0.1*ones(1,10);

%-------------------------------------------------------------------%
%--------------------------- Problem Setup -------------------------%
%-------------------------------------------------------------------%
setup.name                        = 'Tumor-Anti-Angiogenesis-Problem';
setup.functions.continuous        = @tumorAntiAngiogenesisContinuous;
setup.functions.endpoint          = @tumorAntiAngiogenesisEndpoint;
setup.displaylevel                = 2;
setup.auxdata                     = auxdata;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.mesh                        = mesh;
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.scales.method               = 'none';
setup.method                      = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
