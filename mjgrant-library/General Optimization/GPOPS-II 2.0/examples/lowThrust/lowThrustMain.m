% ------------- Low-Thrust Orbit Transfer Example --------------%
% This example is taken verbatim from the following reference:  %
% Betts, J. T., Practical Methods for Optimal Control Using     %
% Nonlinear Programming, SIAM Press, Philadelphia, 2009.        %
% --------------------------------------------------------------%
clc
clear all
close all

% ------------------------------------------------------------- %
%                  Constants and initial conditions             %
% ------------------------------------------------------------- %
T = 4.446618e-3; % [lb]
Isp = 450; % [s]
mu = 1.407645794e16; % [ft^3/s^2]
gs = 32.174; % [ft/s^2]
Re = 20925662.73; % [ft]
J2 = 1082.639e-6;
J3 = -2.565e-6;
J4 = -1.608e-6;

p0 = 21837080.052835; % [ft]
f0 = 0;
g0 = 0;
h0 = -0.25396764647494;
k0 = 0;
L0 = pi;
w0 = 1;

tau = -25;

t0 = 0;
tf = 90000;

% ------------------------------------------------------------- %
%               Set up auxiliary data for problem               %
% ------------------------------------------------------------- %
auxdata.Isp = Isp; % [s]
auxdata.mu = mu; % [ft^3/s^2]
auxdata.gs = gs; % [ft/s^2]
auxdata.T = T; % [lb]
auxdata.Re = Re; % [ft]
auxdata.J2 = J2;
auxdata.J3 = J3;
auxdata.J4 = J4;

% ------------------------------------------------------------- %
%               Generate initial guess for problem              %
% ------------------------------------------------------------- %
initial.p0 = p0;
initial.f0 = f0;
initial.g0 = g0;
initial.h0 = h0;
initial.k0 = k0;
initial.L0 = L0;
initial.w0 = w0;
initial.t0 = t0;

guess.tau = tau;
guess.tf  = tf;

initialguess = lowThrustPropagate(auxdata,initial,guess);

% --------------------------------------------------------------%
%           Set up bounds on state, control, and time           %
% --------------------------------------------------------------%
t0 = 0; tmin = t0; tfmin = 50000; tmax = 100000;
p0 = 21837080.052835; pf = 40007346.015232;
f0 = 0;
g0 = 0;
h0 = -0.25396764647494;
k0 = 0;
L0 = pi;
w0 = 1;
pmin = 20000000; pmax = 60000000;
fmin = -1; fmax = +1;
gmin = -1; gmax = +1;
hmin = -1; hmax = +1;
kmin = -1; kmax = +1;
Lmin = L0; Lmax = 9*2*pi;
wmin = 0.1; wmax = w0;
urmin = -1; urmax = +1;
utmin = -1; utmax = +1;
uhmin = -1; uhmax = +1;
taumin = -50; taumax = 0;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfmin;
bounds.phase.finaltime.upper = tmax;

bounds.phase.initialstate.lower = [p0,f0,g0,h0,k0,L0,w0];
bounds.phase.initialstate.upper = [p0,f0,g0,h0,k0,L0,w0];
bounds.phase.state.lower = [pmin,fmin,gmin,hmin,kmin,Lmin,wmin];
bounds.phase.state.upper = [pmax,fmax,gmax,hmax,kmax,Lmax,wmax];
bounds.phase.finalstate.lower = [pf,fmin,gmin,hmin,kmin,Lmin,wmin];
bounds.phase.finalstate.upper = [pf,fmax,gmax,hmax,kmax,Lmax,wmax];

bounds.phase.control.lower = [urmin,utmin,uhmin];
bounds.phase.control.upper = [urmax,utmax,uhmax];

bounds.parameter.lower = taumin;
bounds.parameter.upper = taumax;

bounds.phase.path.lower = 0;
bounds.phase.path.upper = 0;

bounds.eventgroup.lower = [0.73550320568829^2,0.61761258786099^2,0,-3];
bounds.eventgroup.upper = [0.73550320568829^2,0.61761258786099^2,0,0];

% ----------------------------------------------------- %
% Generate an Initial Guess by Propagating the Dynamics %
% ----------------------------------------------------- %
guess = initialguess;

% ----------------------------------------------------------------- %
% Set Up an Initial Mesh Using the Result from the Propagated Guess %
% ----------------------------------------------------------------- %
mesh.method = 'hp-LiuRao-Legendre';
mesh.tolerance = 1e-5;
mesh.maxiterations = 10;
mesh.colpointsmax = 6;
mesh.colpointsmin = 4;
N = length(guess.fraction);
mesh.phase.colpoints = 4*ones(1,N);
mesh.phase.fraction = guess.fraction;

% ----------------------------------------------------------------------- %
% Set up solver
% ----------------------------------------------------------------------- %
setup.name = 'Betts-Low-Thrust-Orbit-Transfer';
setup.functions.continuous = @lowThrustContinuous;
setup.functions.endpoint = @lowThrustEndpoint;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'adigator';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparseNaN';
setup.derivatives.stepsize = 1e-6;
setup.scales.method = 'automatic-bounds';
setup.mesh = mesh;
setup.method = 'RPM-Differentiation';

% ----------------------------------------------------------------------- %
% Solve problem and extract solution
% ----------------------------------------------------------------------- %
tic;
output = gpops2(setup);
gpopsCPU = toc;

solution = output.result.solution;
