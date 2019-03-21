%--------------- Multiple-Stage Launch Vehicle Ascent Example ------------%
% This example can be found in the following reference:                   %
% Rao, A. V., Benson, D. A., Darby, C. L., Patterson, M. A., Francolin, C.% 
% Sanders, I., and Huntington, G. T., "Algorithm 902: GPOPS, A MATLAB     %
% Software for Solving Multiple-Phase Optimal Control Problems Using the  %
% Gauss Pseudospectral Method," ACM Transactions on Mathematical Software,%
% Vol. 37, No. 2, April-June 2010, Article No. 22, pp. 1-39.              %
% ------------------------------------------------------------------------%
clear all; clc

%-------------------------------------------------------------------------%
%--------------- Provide All Physical Data for Problem -------------------%
%-------------------------------------------------------------------------%
earthRadius         = 6378145;
gravParam           = 3.986012e14;
initialMass         = 301454;
earthRotRate        = 7.29211585e-5;
seaLevelDensity     = 1.225;
densityScaleHeight  = 7200;
g0                  = 9.80665;

scales.length       = 1;
scales.speed        = 1;
scales.time         = 1;
scales.acceleration = 1;
scales.mass         = 1;
scales.force        = 1;
scales.area         = 1;
scales.volume       = 1;
scales.density      = 1;
scales.gravparam    = 1;

if 1,
scales.length       = earthRadius;
scales.speed        = sqrt(gravParam/scales.length);
scales.time         = scales.length/scales.speed;
scales.acceleration = scales.speed/scales.time;
scales.mass         = initialMass;
scales.force        = scales.mass*scales.acceleration;
scales.area         = scales.length^2;
scales.volume       = scales.area.*scales.length;
scales.density      = scales.mass/scales.volume;
scales.pressure     = scales.force/scales.area;
scales.gravparam    = scales.acceleration*scales.length^2;
end

omega               = earthRotRate*scales.time;
auxdata.omegaMatrix = omega*[0 -1 0;1 0 0;0 0 0];
auxdata.mu          = gravParam/scales.gravparam;
auxdata.cd          = 0.5;
auxdata.sa          = 4*pi/scales.area;
auxdata.rho0        = seaLevelDensity/scales.density;
auxdata.H           = densityScaleHeight/scales.length;
auxdata.Re          = earthRadius/scales.length;
auxdata.g0          = g0/scales.acceleration;

lat0                = 28.5*pi/180;         
x0                  = auxdata.Re*cos(lat0);
y0                  = 0;                   
z0                  = auxdata.Re*sin(lat0);
r0                  = [x0 y0 z0];
v0                  = r0*auxdata.omegaMatrix.';
unitr0              = r0/norm(r0,2);
speedrel0           = 5/scales.speed;
v0                  = v0 + speedrel0*unitr0;

btSrb = 75.2/scales.time;
btFirst = 261/scales.time;
btSecond = 700/scales.time;

t0 = 0/scales.time;
t1 = 75.2/scales.time;
t2 = 150.4/scales.time;
t3 = 261/scales.time;
t4 = 961/scales.time;

mTotSrb      = 19290/scales.mass;
mPropSrb     = 17010/scales.mass;
mDrySrb      = mTotSrb-mPropSrb;
mTotFirst    = 104380/scales.mass;
mPropFirst   = 95550/scales.mass;
mDryFirst    = mTotFirst-mPropFirst;
mTotSecond   = 19300/scales.mass;
mPropSecond  = 16820/scales.mass;
mDrySecond   = mTotSecond-mPropSecond;
mPayload     = 4164/scales.mass;
thrustSrb    = 628500/scales.force;
thrustFirst  = 1083100/scales.force;
thrustSecond = 110094/scales.force;
mdotSrb      = mPropSrb/btSrb;
ispSrb       = thrustSrb/(auxdata.g0*mdotSrb);
mdotFirst    = mPropFirst/btFirst;
ispFirst     = thrustFirst/(auxdata.g0*mdotFirst);
mdotSecond   = mPropSecond/btSecond;
ispSecond    = thrustSecond/(auxdata.g0*mdotSecond);

af = 24361140/scales.length;
ef = 0.7308;
incf = 28.5*pi/180;
Omf = 269.8*pi/180;
omf = 130.5*pi/180;
nuguess = 0;
cosincf = cos(incf);
cosOmf = cos(Omf);
cosomf = cos(omf);
oe = [af ef incf Omf omf nuguess];
[rout,vout] = launchoe2rv(oe,auxdata.mu);
rout = rout';
vout = vout';

m10 = mPayload+mTotSecond+mTotFirst+9*mTotSrb;
m1f = m10-(6*mdotSrb+mdotFirst)*t1;
m20 = m1f-6*mDrySrb;
m2f = m20-(3*mdotSrb+mdotFirst)*(t2-t1);
m30 = m2f-3*mDrySrb;
m3f = m30-mdotFirst*(t3-t2);
m40 = m3f-mDryFirst;
m4f = mPayload;

auxdata.thrustSrb    = thrustSrb;
auxdata.thrustFirst  = thrustFirst;
auxdata.thrustSecond = thrustSecond;
auxdata.ispSrb       = ispSrb;
auxdata.ispFirst     = ispFirst;
auxdata.ispSecond    = ispSecond;

rmin = -2*auxdata.Re;
rmax = -rmin;
vmin = -10000/scales.speed;
vmax = -vmin;


%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%
iphase = 1;
bounds.phase(iphase).initialtime.lower = [t0]; 
bounds.phase(iphase).initialtime.upper = [t0]; 
bounds.phase(iphase).finaltime.lower = [t1]; 
bounds.phase(iphase).finaltime.upper = [t1]; 
bounds.phase(iphase).initialstate.lower = [r0(1:3),v0(1:3),m10];    
bounds.phase(iphase).initialstate.upper = [r0(1:3),v0(1:3),m10];    
bounds.phase(iphase).state.lower = [rmin*ones(1,3),vmin*ones(1,3),m1f];
bounds.phase(iphase).state.upper = [rmax*ones(1,3),vmax*ones(1,3),m10];
bounds.phase(iphase).finalstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m1f]; 
bounds.phase(iphase).finalstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m10]; 
bounds.phase(iphase).control.lower = -10*ones(1,3);
bounds.phase(iphase).control.upper = +10*ones(1,3);
bounds.phase(iphase).path.lower  = [1];
bounds.phase(iphase).path.upper  = [1];
guess.phase(iphase).time = [t0; t1];
guess.phase(iphase).state(:,1) = [r0(1); r0(1)];
guess.phase(iphase).state(:,2) = [r0(2); r0(2)];
guess.phase(iphase).state(:,3) = [r0(3); r0(3)];
guess.phase(iphase).state(:,4) = [v0(1); v0(1)];
guess.phase(iphase).state(:,5) = [v0(2); v0(2)];
guess.phase(iphase).state(:,6) = [v0(3); v0(3)];
guess.phase(iphase).state(:,7) = [m10; m1f];
guess.phase(iphase).control(:,1) = [0; 0];
guess.phase(iphase).control(:,2) = [1; 1];
guess.phase(iphase).control(:,3) = [0; 0];

iphase = 2;
bounds.phase(iphase).initialtime.lower = [t1]; 
bounds.phase(iphase).initialtime.upper = [t1]; 
bounds.phase(iphase).finaltime.lower = [t2]; 
bounds.phase(iphase).finaltime.upper = [t2]; 
bounds.phase(iphase).initialstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m2f];
bounds.phase(iphase).initialstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m20];
bounds.phase(iphase).state.lower = [rmin*ones(1,3),vmin*ones(1,3),m2f];
bounds.phase(iphase).state.upper = [rmax*ones(1,3),vmax*ones(1,3),m20];
bounds.phase(iphase).finalstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m2f]; 
bounds.phase(iphase).finalstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m20]; 
bounds.phase(iphase).control.lower = -10*ones(1,3);
bounds.phase(iphase).control.upper = +10*ones(1,3);
bounds.phase(iphase).path.lower  = 1;
bounds.phase(iphase).path.upper  = 1;
guess.phase(iphase).time = [t1; t2];
guess.phase(iphase).state(:,1) = [r0(1); r0(1)];
guess.phase(iphase).state(:,2) = [r0(2); r0(2)];
guess.phase(iphase).state(:,3) = [r0(3); r0(3)];
guess.phase(iphase).state(:,4) = [v0(1); v0(1)];
guess.phase(iphase).state(:,5) = [v0(2); v0(2)];
guess.phase(iphase).state(:,6) = [v0(3); v0(3)];
guess.phase(iphase).state(:,7) = [m10; m1f];
guess.phase(iphase).control(:,1) = [0; 0];
guess.phase(iphase).control(:,2) = [1; 1];
guess.phase(iphase).control(:,3) = [0; 0];

iphase = 3;
bounds.phase(iphase).initialtime.lower = [t2]; 
bounds.phase(iphase).initialtime.upper = [t2]; 
bounds.phase(iphase).finaltime.lower = [t3]; 
bounds.phase(iphase).finaltime.upper = [t3]; 
bounds.phase(iphase).initialstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m3f];
bounds.phase(iphase).initialstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m30];
bounds.phase(iphase).state.lower = [rmin*ones(1,3),vmin*ones(1,3),m3f];
bounds.phase(iphase).state.upper = [rmax*ones(1,3),vmax*ones(1,3),m30];
bounds.phase(iphase).finalstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m3f]; 
bounds.phase(iphase).finalstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m30]; 
bounds.phase(iphase).control.lower = -10*ones(1,3);
bounds.phase(iphase).control.upper = +10*ones(1,3);
bounds.phase(iphase).path.lower  = 1;
bounds.phase(iphase).path.upper  = 1;
guess.phase(iphase).time       = [t2; t3];
guess.phase(iphase).state(:,1) = [rout(1); rout(1)];
guess.phase(iphase).state(:,2) = [rout(2); rout(2)];
guess.phase(iphase).state(:,3) = [rout(3); rout(3)];
guess.phase(iphase).state(:,4) = [vout(1); vout(1)];
guess.phase(iphase).state(:,5) = [vout(2); vout(2)];
guess.phase(iphase).state(:,6) = [vout(3); vout(3)];
guess.phase(iphase).state(:,7) = [m30; m3f];
guess.phase(iphase).control(:,1) = [1; 1];
guess.phase(iphase).control(:,2) = [0; 0];
guess.phase(iphase).control(:,3) = [0; 0];

iphase = 4;
bounds.phase(iphase).initialtime.lower = [t3]; 
bounds.phase(iphase).initialtime.upper = [t3]; 
bounds.phase(iphase).finaltime.lower = [t3]; 
bounds.phase(iphase).finaltime.upper = [t4]; 
bounds.phase(iphase).initialstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m4f];
bounds.phase(iphase).initialstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m40];
bounds.phase(iphase).state.lower = [rmin*ones(1,3),vmin*ones(1,3),m4f];
bounds.phase(iphase).state.upper = [rmax*ones(1,3),vmax*ones(1,3),m40];
bounds.phase(iphase).finalstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m4f]; 
bounds.phase(iphase).finalstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m40]; 
bounds.phase(iphase).control.lower = -10*ones(1,3);
bounds.phase(iphase).control.upper = +10*ones(1,3);
bounds.phase(iphase).path.lower  = 1;
bounds.phase(iphase).path.upper  = 1;
guess.phase(iphase).time    = [t3; t4];
guess.phase(iphase).state(:,1) = [rout(1) rout(1)];
guess.phase(iphase).state(:,2) = [rout(2) rout(2)];
guess.phase(iphase).state(:,3) = [rout(3) rout(3)];
guess.phase(iphase).state(:,4) = [vout(1) vout(1)];
guess.phase(iphase).state(:,5) = [vout(2) vout(2)];
guess.phase(iphase).state(:,6) = [vout(3) vout(3)];
guess.phase(iphase).state(:,7) = [m40; m4f];
guess.phase(iphase).control(:,1) = [1; 1];
guess.phase(iphase).control(:,2) = [0; 0];
guess.phase(iphase).control(:,3) = [0; 0];

%-------------------------------------------------------------------------%
%------------- Set up Event Constraints That Link Phases -----------------%
%-------------------------------------------------------------------------%
bounds.eventgroup(1).lower = [zeros(1,6), -6*mDrySrb, 0];
bounds.eventgroup(1).upper = [zeros(1,6), -6*mDrySrb, 0];
bounds.eventgroup(2).lower = [zeros(1,6), -3*mDrySrb, 0];
bounds.eventgroup(2).upper = [zeros(1,6), -3*mDrySrb, 0];
bounds.eventgroup(3).lower = [zeros(1,6), -mDryFirst, 0];
bounds.eventgroup(3).upper = [zeros(1,6), -mDryFirst, 0];

%-------------------------------------------------------------------------%
%----------- Set up Event Constraints That Define Final Orbit ------------%
%-------------------------------------------------------------------------%
bounds.eventgroup(4).lower = [af, ef, incf, Omf, omf];
bounds.eventgroup(4).upper = [af, ef, incf, Omf, omf];

%-------------------------------------------------------------------------%
%-------------- Provide an Initial Mesh in Each Phase --------------------%
%-------------------------------------------------------------------------%
for i=1:4
  meshphase(i).colpoints = 4*ones(1,10);
  meshphase(i).fraction = 0.1*ones(1,10);
end

%-------------------------------------------------------------------------%
%----------- Assemble All Information into Setup Structure ---------------%
%-------------------------------------------------------------------------%
setup.name = 'Launch-Vehicle-Ascent-Problem';
setup.functions.continuous = @launchContinuous;
setup.functions.endpoint = @launchEndpoint;
setup.mesh.phase = meshphase;
setup.nlp.solver = 'ipopt';
setup.nlp.snoptoptions.tolerance = 1e-7;
setup.bounds = bounds;
setup.guess = guess;
setup.auxdata = auxdata;
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparseNaN';
% setup.scales.method = 'automatic-bounds';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-6;
setup.method = 'RPM-Differentiation';

%-------------------------------------------------------------------------%
%---------------------- Solve Problem using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
