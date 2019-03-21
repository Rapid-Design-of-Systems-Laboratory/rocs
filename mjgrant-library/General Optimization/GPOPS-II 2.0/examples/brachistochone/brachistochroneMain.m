%---------------------------------------------------%
% Classical Brachistochrone Problem:                %
%---------------------------------------------------%
% The problem solved here is given as follows:      %
%   Minimize t_f                                    %
% subject to the dynamic constraints                %
%    dx/dt = v*sin(u)                               %
%    dy/dt = v*cos(u)                               %
%    dv/dt = g*cos(u)                               %
% and the boundary conditions                       %
%    x(0) = 0, y(0) = 0, v(0) = 0                   %
%    x(t_f) = 2, y(t_f) = 2, v(t_f) = FREE          %
%---------------------------------------------------%

clear all; close all; clc

auxdata.g = 10;
t0 = 0; 
tfmin = 0; tfmax = 10;
x0 = 0; y0 = 0; v0 = 0;
xf = 2; yf = 2;
xmin = 0; xmax = 10;
ymin = 0; ymax = 10;
vmin = -50; vmax = 50;
umin = -pi/2; umax = pi/2;

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfmin; 
bounds.phase.finaltime.upper = tfmax;
bounds.phase.initialstate.lower = [x0,y0,v0]; 
bounds.phase.initialstate.upper = [x0,y0,v0]; 
bounds.phase.state.lower = [xmin,ymin,vmin]; 
bounds.phase.state.upper = [xmax,ymax,vmax]; 
bounds.phase.finalstate.lower = [xf,yf,vmin]; 
bounds.phase.finalstate.upper = [xf,yf,vmax]; 
bounds.phase.control.lower = umin; 
bounds.phase.control.upper = umax;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [t0; tfmax]; 
guess.phase.state   = [[x0; xf],[y0; yf],[v0; v0]];
guess.phase.control = [0; pi/2];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-PattersonRao';
mesh.tolerance    = 1e-6;
mesh.maxiterations = 45;
mesh.colpointsmin = 4;
mesh.colpointsmax = 10;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                        = 'Brachistochrone-Problem';
setup.functions.continuous        = @brachistochroneContinuous;
setup.functions.endpoint          = @brachistochroneEndpoint;
setup.auxdata                     = auxdata;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.mesh                        = mesh; 
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.method                      = 'RPM-Differentiation';

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOP2 ---------------------%
%-------------------------------------------------------------------------%
output   = gpops2(setup);
solution = output.result.solution;
