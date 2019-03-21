%--------------- Multiple-Stage Launch Vehicle Ascent Example ------------%
% This example can be found in the following reference:                   %
% Rao, A. V., Benson, D. A., Darby, C. L., Patterson, M. A., Francolin,   % 
% C., Sanders, I., and Huntington, G. T., “Algorithm 902: GPOPS, A MATLAB %
% Software for Solving Multiple-Phase Optimal Control Problems Using The  %
% Gauss Pseudospectral Method,” ACM Transactions on Mathematical Software,%
% Vol. 37, No. 2, April-June 2010, Article No. 22, pp. 1-39.              %
% ------------------------------------------------------------------------%
clear all; clc

h0 = 0;
v0 = 0;
m0 = 3;
mf = 1;
hmin = 0;
hmax = 30000;
vmin = -1500;
vmax = 1500;
mmin = 0.2*m0;
mmax = m0;
t0 = 0;
tfMin = 0;
tfMax = 500;

auxdata.g0 = 32.174;
auxdata.rho0 = 0.002378;
auxdata.H = 23800;
auxdata.csqrd = 3.264*auxdata.g0*auxdata.H;
auxdata.c = sqrt(auxdata.csqrd);
Tmax = 2*m0*auxdata.g0;
auxdata.dragk = 0.7110*Tmax/auxdata.csqrd;

%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%
iphase = 1;
bounds.phase(iphase).initialtime.lower = t0;
bounds.phase(iphase).initialtime.upper = t0;
bounds.phase(iphase).finaltime.lower = tfMin;
bounds.phase(iphase).finaltime.upper = tfMax;
bounds.phase(iphase).initialstate.lower = [h0, v0, m0];
bounds.phase(iphase).initialstate.upper = [h0, v0, m0];
bounds.phase(iphase).state.lower = [hmin, vmin, mmin];
bounds.phase(iphase).state.upper = [hmax, vmax, mmax];
bounds.phase(iphase).finalstate.lower = [hmin, vmin, mmin];
bounds.phase(iphase).finalstate.upper = [hmax, vmax, mmax];
bounds.phase(iphase).control.lower = Tmax;
bounds.phase(iphase).control.upper = Tmax;
guess.phase(iphase).time =  [t0; tfMax];
guess.phase(iphase).state(:,1) = [h0; h0];
guess.phase(iphase).state(:,2) = [v0; v0];
guess.phase(iphase).state(:,3) = [m0; mf];
guess.phase(iphase).control = [0; Tmax];

iphase = 2;
bounds.phase(iphase).initialtime.lower = tfMin;
bounds.phase(iphase).initialtime.upper = tfMax;
bounds.phase(iphase).finaltime.lower = tfMin;
bounds.phase(iphase).finaltime.upper = tfMax;
bounds.phase(iphase).initialstate.lower = [hmin, vmin, mmin];
bounds.phase(iphase).initialstate.upper = [hmax, vmax, mmax];
bounds.phase(iphase).state.lower = [hmin, vmin, mmin];
bounds.phase(iphase).state.upper = [hmax, vmax, mmax];
bounds.phase(iphase).finalstate.lower = [hmin, vmin, mmin];
bounds.phase(iphase).finalstate.upper = [hmax, vmax, mmax];
bounds.phase(iphase).control.lower = 0;
bounds.phase(iphase).control.upper = Tmax;
bounds.phase(iphase).path.lower = 0;
bounds.phase(iphase).path.upper = 0;
guess.phase(iphase).time =  [tfMin; tfMax];
guess.phase(iphase).state(:,1) = [h0; hmax];
guess.phase(iphase).state(:,2) = [vmax; vmax];
guess.phase(iphase).state(:,3) = [m0; mf];
guess.phase(iphase).control(:,1) = [0; Tmax];

iphase = 3;
bounds.phase(iphase).initialtime.lower = tfMin;
bounds.phase(iphase).initialtime.upper = tfMax;
bounds.phase(iphase).finaltime.lower = tfMin;
bounds.phase(iphase).finaltime.upper = tfMax;
bounds.phase(iphase).initialstate.lower = [hmin, vmin, mmin];
bounds.phase(iphase).initialstate.upper = [hmax, vmax, mmax];
bounds.phase(iphase).state.lower = [hmin, vmin, mmin];
bounds.phase(iphase).state.upper = [hmax, vmax, mmax];
bounds.phase(iphase).finalstate.lower = [hmin, vmin, mf];
bounds.phase(iphase).finalstate.upper = [hmax, vmax, mf];
bounds.phase(iphase).control.lower = 0;
bounds.phase(iphase).control.upper = 0;
guess.phase(iphase).time =  [tfMin; tfMax];
guess.phase(iphase).state(:,1) = [hmin; hmax];
guess.phase(iphase).state(:,2) = [vmax; vmax];
guess.phase(iphase).state(:,3) = [m0; mf];
guess.phase(iphase).control(:,1) = [0; Tmax];

%-------------------------------------------------------------------------%
%------------- Set up Event Constraints That Link Phases -----------------%
%-------------------------------------------------------------------------%
bounds.eventgroup(1).lower = zeros(1,4);
bounds.eventgroup(1).upper = zeros(1,4);
bounds.eventgroup(2).lower = 0;
bounds.eventgroup(2).upper = 0;
bounds.eventgroup(3).lower = zeros(1,4);
bounds.eventgroup(3).upper = zeros(1,4);
bounds.eventgroup(4).lower = 0.1*ones(1,3);
bounds.eventgroup(4).upper = 1000*ones(1,3);

%-------------------------------------------------------------------------%
%----------- Assemble All Information into Setup Structure ---------------%
%-------------------------------------------------------------------------%
setup.name = 'Goddard-Rocket-Problem';
setup.functions.continuous = @goddardRocketContinuous;
setup.functions.endpoint = @goddardRocketEndpoint;
setup.nlp.solver = 'ipopt';
setup.bounds = bounds;
setup.guess = guess;
setup.auxdata = auxdata;
% setup.derivatives.supplier = 'adigator';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparseNaN';
setup.scales.method = 'automatic-bounds';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-6;
setup.method = 'RPM-Differentiation';
setup.derivatives.supplier = 'adigator';
%-------------------------------------------------------------------------%
%---------------------- Solve Problem using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
