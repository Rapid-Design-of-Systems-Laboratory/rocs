%------------------- Two-Strain Tuberculosis Model ------------------%
% This problem is taken from the following reference:                %
% Betts, J. T., "Practical Methods for Optimal Control and           %
% Estimation Using Nonlinear Programming," SIAM Press, Philadelphia, %
% PA, 2009.                                                          %
%--------------------------------------------------------------------%
clear all
close all
clc

%-------------------------------------------------------------------%
%-------------- Provide Auxiliary Data for Problem -----------------%
%-------------------------------------------------------------------%
auxdata.beta1 = 13;
auxdata.beta2 = 13;
auxdata.mu = 0.0143;
auxdata.d1 = 0;
auxdata.d2 = 0;
auxdata.k1 = 0.5;
auxdata.k2 = 1;
auxdata.r1 = 2;
auxdata.r2 = 1;
auxdata.p = 0.4;
auxdata.q = 0.1;
auxdata.Npop = 30000;
auxdata.betas = 0.029;
auxdata.B1 = 50;
auxdata.B2 = 500;
auxdata.lam = auxdata.mu*auxdata.Npop;
auxdata.m0 = 1;
auxdata.dm = 0.0749;

%-------------------------------------------------------------------%
%----------- Set up Bounds for Optimal Control Problem -------------%
%-------------------------------------------------------------------%
t0 = 0;
tf = 5;
S0 = 76*auxdata.Npop/120;
T0 = auxdata.Npop/120;
L10 = 36*auxdata.Npop/120;
L20 = 2*auxdata.Npop/120;
I10 = 4*auxdata.Npop/120;
I20 = 1*auxdata.Npop/120;
Nmin = 0;
Nmax = 30000;
u1min = 0.05;
u1max = 0.95;
u2min = 0.05;
u2max = 0.95;

bounds.phase.initialstate.lower = [S0, T0, L10, L20, I10, I20];
bounds.phase.initialstate.upper = [S0, T0, L10, L20, I10, I20];
bounds.phase.state.lower = [Nmin, Nmin, Nmin, Nmin, Nmin, Nmin];
bounds.phase.state.upper = [Nmax, Nmax, Nmax, Nmax, Nmax, Nmax];
bounds.phase.finalstate.lower = [Nmin, Nmin, Nmin, Nmin, Nmin, Nmin];
bounds.phase.finalstate.upper = [Nmax, Nmax, Nmax, Nmax, Nmax, Nmax];
bounds.phase.control.lower = [u1min, u2min];
bounds.phase.control.upper = [u1max, u2max];
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.integral.lower = [0];
bounds.phase.integral.upper = [10000];
bounds.phase.path.lower  = [0];
bounds.phase.path.upper  = [0];

%-------------------------------------------------------------------%
%------------ Provide an Initial Guess of the Solution -------------%
%-------------------------------------------------------------------%
timeGuess = [t0; tf];
SGuess    = [S0; S0];
TGuess    = [T0; S0];
L1Guess   = [L10; L10];
L2Guess   = [L20; L20];
I1Guess   = [I10; I10];
I2Guess   = [I20; I20];
u1Guess   = [0.95; 0.95];
u2Guess   = [0.95; 0.95];
guess.phase.time    = [timeGuess];
guess.phase.state   = [SGuess, TGuess, L1Guess, L2Guess, I1Guess, I2Guess];
guess.phase.control = [u1Guess, u2Guess];
guess.phase.integral = 6000;

%-------------------------------------------------------------------%
%------------ Provide an Initial Mesh for the Solution -------------%
%-------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-6;
mesh.maxiterations   = 10;
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 10;
N = 10;
mesh.phase.colpoints = 4*ones(1,N);
mesh.phase.fraction   = ones(1,N)/N;

setup.name                        = 'Tuberculosis-Optimal-Control-Problem';
setup.functions.continuous        = @tuberculosisContinuous;
setup.functions.endpoint          = @tuberculosisEndpoint;
setup.displaylevel                = 2;
setup.auxdata                     = auxdata;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.mesh                        = mesh;
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'adigator';
setup.derivatives.derivativelevel = 'second';
setup.scales.method               = 'automatic-bounds';
setup.method                      = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%--------- Solve Problem with GPOPS2 and Extract Solution ----------%
%-------------------------------------------------------------------%
output = gpops2(setup);
