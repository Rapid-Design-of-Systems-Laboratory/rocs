%------------------- Dynamic Soaring Problem -----------------------------%
% This example is taken from the following reference:                     %
% Zhao, Y. J., "Optimal Pattern of Glider Dynamic Soaring," Optimal       %
% Control Applications and Methods, Vol. 25, 2004, pp. 67-89.             %
%-------------------------------------------------------------------------%
clear all
clc

%-------------------------------------------------------------------------%
%------------------ Provide Auxiliary Data for Problem -------------------%
%-------------------------------------------------------------------------%
auxdata.rho  = 0.002378;
auxdata.CD0  = 0.00873; 
auxdata.K    = 0.045;
auxdata.g    = 32.2;
auxdata.m    = 5.6;
auxdata.S    = 45.09703;
auxdata.mu   = 3.986e14;
auxdata.mgos = auxdata.m*auxdata.g/auxdata.S;
auxdata.Emax = (1/(4*auxdata.K*auxdata.CD0))^0.5;
auxdata.W0   = 0;
auxdata.lmin = -2;
auxdata.lmax = 5;

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
t0 = 0; x0 = 0; y0 = 0; z0 = 0; v0 = 100;

%-------------------------------------------------------------------%
%----------------------- Limits on Variables -----------------------%
%-------------------------------------------------------------------%
tfmin    = 1;          tfmax    = 30;
xmin     = -1000;      xmax     = +1000;
ymin     = -1000;      ymax     = +1000;
zmin     =  0;         zmax     = +1000;
vmin     = +10;        vmax     = +350;
gammamin = -75*pi/180; gammamax =  75*pi/180;
psimin   = -3*pi;      psimax   = +pi/2;
betamin  = 0.005;      betamax  = 0.15;
CLmin    = -0.5;       CLmax    = 1.5;
Phimin   = -75/180*pi; Phimax   =  75/180*pi;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase.initialtime.lower  = t0;
bounds.phase.initialtime.upper  = t0;
bounds.phase.finaltime.lower    = tfmin;
bounds.phase.finaltime.upper    = tfmax;
bounds.phase.initialstate.lower = [x0, y0, z0, vmin, gammamin, psimin];
bounds.phase.initialstate.upper = [x0, y0, z0, vmax, gammamax, psimax];
bounds.phase.state.lower        = [xmin, ymin, zmin, vmin, gammamin, psimin];
bounds.phase.state.upper        = [xmax, ymax, zmax, vmax, gammamax, psimax];
bounds.phase.finalstate.lower   = [x0, y0, z0, vmin, gammamin, psimin];
bounds.phase.finalstate.upper   = [x0, y0, z0, vmax, gammamax, psimax];
bounds.phase.control.lower      = [CLmin, Phimin];
bounds.phase.control.upper      = [CLmax, Phimax];
bounds.phase.path.lower         = auxdata.lmin;
bounds.phase.path.upper         = auxdata.lmax;
bounds.eventgroup(1).lower      = [0, 0, -2*pi];
bounds.eventgroup(1).upper      = [0, 0, -2*pi];
bounds.parameter.lower          = betamin;
bounds.parameter.upper          = betamax;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
N                               = 100;
CL0                             = CLmax;
tGuess                          = linspace(0,24,N).';
xguess                          = 500*cos(2*pi*tGuess/24)-500;
yguess                          = 300*sin(2*pi*tGuess/24);
zguess                          = -400*cos(2*pi*tGuess/24)+400;
vguess                          = 0.8*v0*(1.5+cos(2*pi*tGuess/24));
gammaguess                      = pi/6*sin(2*pi*tGuess/24);
psiguess                        = -1-tGuess/4;
CLguess                         = CL0*ones(N,1)/3;
phiguess                        = -ones(N,1);
betaguess                       = 0.08;
guess.phase.time                = tGuess;
guess.phase.state               = [xguess, yguess, zguess, vguess, gammaguess, psiguess];
guess.phase.control             = [CLguess, phiguess];
guess.parameter                 = [betaguess];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.maxiterations              = 10;
mesh.method                     = 'hp-LiuRao';
mesh.tolerance                  = 1e-6;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                             = 'Dynamic-Soaring-Problem';
setup.functions.continuous             = @dynamicSoaringContinuous;
setup.functions.endpoint               = @dynamicSoaringEndpoint;
setup.nlp.solver                       = 'ipopt';
setup.nlp.ipoptoptions.linear_solver   = 'ma57';
setup.displaylevel                     = 2;
setup.auxdata                          = auxdata;
setup.bounds                           = bounds;
setup.guess                            = guess;
setup.mesh                             = mesh;
setup.derivatives.supplier             = 'adigator';
setup.derivatives.derivativelevel      = 'second';
setup.scales.method                    = 'automatic-bounds';
setup.method                           = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------ Solve Problem Using GPOPS-II--------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
solution = output.result.solution;
