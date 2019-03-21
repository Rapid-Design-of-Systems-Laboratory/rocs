%-------------------- Kinetic Batch Reactor Problem ----------------------%
% This problem is taken from the following reference:                     %
% Betts, J. T., Practical Methods for Optimal Control and Estimation      %
% Using Nonlinear Programming, SIAM Press, Philadelphia, PA, 2009.        %
%-------------------------------------------------------------------------%
clear all
close all
clc

%-------------------------------------%
%            Problem Setup            %
%-------------------------------------%
chr2min  = 60;               % Conversion from Hours to Minutes
cmin2hr  = 1/chr2min;        % Conversion from Minutes to Hours
cmin2sec = 60;               % Conversion from Minutes to Seconds
csec2min = 1/cmin2sec;       % Conversion from Seconds to Minutes
chr2sec  = chr2min*cmin2sec; % Conversion from Hours to Seconds
csec2hr  = 1/chr2sec;        % Conversion from Seconds to Hours
csec2ms  = 1000;             % Conversion from Seconds to Milliseconds
cmin2ms  = cmin2sec*csec2ms; % Conversion from Minutes to Milliseconds
chr2ms   = chr2min*cmin2ms;  % Conversion from Hours to Milliseconds
auxdata.chr2min = chr2min;
auxdata.cmin2hr = cmin2hr;
auxdata.cmin2sec = cmin2sec;
auxdata.csec2min = csec2min;
auxdata.chr2sec = chr2sec;
auxdata.csec2hr = csec2hr;
auxdata.csec2ms = csec2ms;
auxdata.cmin2ms  = cmin2ms;
auxdata.chr2ms  = chr2ms;

k1hat      = 1.3708e12;  % kg/gmol/hr
kminus1hat = 1.6215e20;  % kg/gmol/hr
k2hat      = 5.2282e12;  % kg/gmol/hr
beta1      = 9.2984e3;   % K
betaminus1 = 1.3108e4;   % K
beta2      = 9.5999e3;   % K
K1         = 2.575e-16;  % gmol/kg
K2         = 4.876e-14;  % gmol/kg
K3         = 1.7884e-16; % gmol/kg
a          = 2;          % gmol/kg/hr^2

auxdata.phase(1).k1hat      = k1hat;
auxdata.phase(1).kminus1hat = kminus1hat;
auxdata.phase(1).k2hat      = k2hat;
auxdata.phase(1).a          = a;

auxdata.phase(2).k1hat      = k1hat;
auxdata.phase(2).kminus1hat = kminus1hat;
auxdata.phase(2).k2hat      = k2hat;
auxdata.phase(2).a          = a;

auxdata.phase(3).k1hat      = k1hat;
auxdata.phase(3).kminus1hat = kminus1hat;
auxdata.phase(3).k2hat      = k2hat;
auxdata.phase(3).a          = a;

auxdata.beta1               = beta1;
auxdata.betaminus1          = betaminus1;
auxdata.beta2               = beta2;
auxdata.K1                  = K1;
auxdata.K2                  = K2;
auxdata.K3                  = K3;

tfupphr   = 10;   % Upper Bound on Final Time (Hours)
tfuppmin  = tfupphr*chr2min;
tfuppms   = tfupphr*chr2ms;
tfuppsec  = tfupphr*chr2sec;
t10hr     = 0;
t10sec    = t10hr*chr2sec;
t10ms     = t10hr*chr2ms;
t1fhr     = 0.01;
t1fsec    = t1fhr*chr2sec;
t1fms     = t1fhr*chr2ms;
t20hr     = t1fhr;
t20lowhr  = t20hr;
t20upphr = tfupphr;
t20lowms  = t20hr*chr2ms;
t20uppms = tfupphr*chr2ms;
t20lowsec = t20hr*chr2sec;
t20uppsec = tfupphr*chr2sec;
t20lowmin = t20hr*chr2min;
t20uppmin = tfupphr*chr2min;
t2flowhr  = t20lowhr;
t2fupphr  = tfupphr;
t2flowms  = t20lowms;
t2fuppms  = tfupphr*chr2ms;
t2flowsec = t20lowsec;
t2fuppsec = tfupphr*chr2sec;
t2flowmin = t20lowmin;
t2fuppmin = tfupphr*chr2min;
t30lowms   = t20lowms;
t30uppms   = tfupphr*chr2ms;
t30lowsec  = t20lowsec;
t30uppsec  = tfupphr*chr2sec;
t30lowmin  = t20lowmin;
t30uppmin  = tfupphr*chr2min;
t3flowms   = t30lowms;
t3fuppms   = tfupphr*chr2ms;
t3flowsec  = t30lowsec;
t3fuppsec  = tfupphr*chr2sec;
t3flowmin  = t30lowmin;
t3fuppmin  = tfupphr*chr2min;
t30lowhr   = t20lowmin/chr2min;
t30upphr  = tfupphr;
t3flowhr  = t30lowhr;
t3fupphr  = tfupphr;

%-------------------------------------%
%          Boundary Conditions        %
%-------------------------------------%
y10 = 1.5776; y20 = 8.32; y30 = 0; y40 = 0; y50 = 0;
y1min = 0; y1max = 2;
y2min = 5; y2max = 10;
y3min = 0; y3max = 1;
y4min = 0; y4max = 1.5;
y5min = 0; y5max = 1.5;
% y6min = 0; y6max = 0.02;
y6min = 0; y6max = 0.1;
% u1min = 6; u1max = 12;
u1min = 0; u1max = 15;
u2min = 0; u2max = 0.02;
u3min = 0; u3max = 2e-5;
u4min = 0; u4max = 2e-5;
u5min = 293.15; u5max = 393.15;
pmin = 0; pmax = 0.0262;

%-----------------------------------------------%
% Bounds on Time, State, and Control in Phase 1 %
%-----------------------------------------------%
iphase = 1;
% Numsec = 1;
% meshphase(iphase).colpoints = 15.*ones(1,Numsec);
% meshphase(iphase).fraction   = 1/Numsec.*ones(1,Numsec);
Numsec = 10;
meshphase(iphase).colpoints = 3.*ones(1,Numsec);
meshphase(iphase).fraction   = 1/Numsec.*ones(1,Numsec);
bounds.phase(iphase).initialtime.lower = t10hr;
bounds.phase(iphase).initialtime.upper = t10hr;
bounds.phase(iphase).finaltime.lower   = t1fhr;
bounds.phase(iphase).finaltime.upper   = t1fhr;
bounds.phase(iphase).initialstate.lower = [y10,y20,y30,y40,y50,y6min];
bounds.phase(iphase).initialstate.upper = [y10,y20,y30,y40,y50,y6max];
bounds.phase(iphase).state.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).state.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).finalstate.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).finalstate.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).control.lower = [u1min,u2min,u3min,u4min,u5min];
bounds.phase(iphase).control.upper = [u1max,u2max,u3max,u4max,u5max];
bounds.phase(iphase).path.lower = [0,0,0,0,0];
bounds.phase(iphase).path.upper = [0,0,0,0,100];

%--------------------------------------%
% Initial Guess of Solution in Phase 1 %
%--------------------------------------%
tguess = [t10hr; t1fhr];
y1guess = [y10; 0.5];
y2guess = [y20; 6];
y3guess = [y30; 0.6];
y4guess = [0; 0.5];
y5guess = [0; 0.5];
y6guess = [0.013; 0.013];
u1guess = [7; 10];
u2guess = [0; 0];
u3guess = [0; 1e-5];
u4guess = [0; 1e-5];
u5guess = [373; 393.15];
guess.phase(iphase).state   = [y1guess,y2guess,y3guess,y4guess,y5guess,y6guess];
guess.phase(iphase).control = [u1guess,u2guess,u3guess,u4guess,u5guess];
guess.phase(iphase).time    = tguess;

%-----------------------------------------------%
% Bounds on Time, State, and Control in Phase 2 %
%-----------------------------------------------%
iphase = 2;
% Numsec = 1;
% meshphase(iphase).colpoints = 15.*ones(1,Numsec);
% meshphase(iphase).fraction   = 1/Numsec.*ones(1,Numsec);
Numsec = 10;
meshphase(iphase).colpoints = 3.*ones(1,Numsec);
meshphase(iphase).fraction   = 1/Numsec.*ones(1,Numsec);
bounds.phase(iphase).initialtime.lower = t20lowhr;
bounds.phase(iphase).initialtime.upper = t20upphr;
bounds.phase(iphase).finaltime.lower   = t2flowhr;
bounds.phase(iphase).finaltime.upper   = t2fupphr;
bounds.phase(iphase).initialstate.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).initialstate.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).state.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).state.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).finalstate.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).finalstate.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).control.lower = [u1min,u2min,u3min,u4min,u5min];
bounds.phase(iphase).control.upper = [u1max,u2max,u3max,u4max,u5max];
bounds.phase(iphase).path.lower = [0,0,0,0,0];
bounds.phase(iphase).path.upper = [0,0,0,0,100];

%--------------------------------------%
% Initial Guess of Solution in Phase 2 %
%--------------------------------------%
tguess = [t20lowhr; tfupphr/4];
y1guess = [0.5; 0.5];
y2guess = [6; 6];
y3guess = [0.6; 0.6];
y4guess = [0.5; 0.5];
y5guess = [0.5; 0.5];
y6guess = [0.013; 0.013];
u1guess = [10; 10];
u2guess = [0; 0];
u3guess = [1e-5; 1e-5];
u4guess = [1e-5; 1e-5];
u5guess = [393.15; 393.15];
guess.phase(iphase).state   = [y1guess,y2guess,y3guess,y4guess,y5guess,y6guess];
guess.phase(iphase).control = [u1guess,u2guess,u3guess,u4guess,u5guess];
guess.phase(iphase).time    = tguess;

%-----------------------------------------------%
% Bounds on Time, State, and Control in Phase 3 %
%-----------------------------------------------%
iphase = 3;
% Numsec = 1;
% meshphase(iphase).colpoints = 15.*ones(1,Numsec);
% meshphase(iphase).fraction   = 1/Numsec.*ones(1,Numsec);
Numsec = 10;
meshphase(iphase).colpoints = 3.*ones(1,Numsec);
meshphase(iphase).fraction   = 1/Numsec.*ones(1,Numsec);
y4fmin = 1;
y4fmax = y4max;
bounds.phase(iphase).initialtime.lower = t30lowhr;
bounds.phase(iphase).initialtime.upper = t30upphr;
bounds.phase(iphase).finaltime.lower   = t3flowhr;
bounds.phase(iphase).finaltime.upper   = t3fupphr;
bounds.phase(iphase).initialstate.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).initialstate.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).state.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).state.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).finalstate.lower = [y1min,y2min,y3min,y4min,y5min,y6min];
bounds.phase(iphase).finalstate.upper = [y1max,y2max,y3max,y4max,y5max,y6max];
bounds.phase(iphase).control.lower = [u1min,u2min,u3min,u4min,u5min];
bounds.phase(iphase).control.upper = [u1max,u2max,u3max,u4max,u5max];
bounds.phase(iphase).path.lower = [0,0,0,0];
bounds.phase(iphase).path.upper = [0,0,0,0];

%--------------------------------------%
% Initial Guess of Solution in Phase 3 %
%--------------------------------------%
tguess = [tfupphr/4; tfupphr];
y1guess = [0.5; 0.1];
y2guess = [6; 5];
y3guess = [6; 6];
y4guess = [0.5; 0.9];
y5guess = [0.5; 0.9];
y6guess = [0.013; 0.013];
u1guess = [10; 10];
u2guess = [0; 0];
u3guess = [1e-5; 1e-5];
u4guess = [1e-5; 1e-5];
u5guess = [393.15; 393.15];
guess.phase(iphase).state   = [y1guess,y2guess,y3guess,y4guess,y5guess,y6guess];
guess.phase(iphase).control = [u1guess,u2guess,u3guess,u4guess,u5guess];
guess.phase(iphase).time    = tguess;

%---------------------%
% Bounds on Parameter %
%---------------------%
bounds.parameter.lower = 0;
bounds.parameter.upper = pmax;

%----------------------------%
% Initial Guess on Parameter %
%----------------------------%
guess.parameter = 0.013;

%----------------------------------------------------------------------------%
% Event Group 1: Linkage Constraint on State and Time Between Phases 1 and 2 %
%----------------------------------------------------------------------------%
bounds.eventgroup(1).lower = [zeros(1,7)];
bounds.eventgroup(1).upper = [zeros(1,7)];
%----------------------------------------------------------------------------%
% Event Group 2: Linkage Constraint on State and Time Between Phases 2 and 3 %
%----------------------------------------------------------------------------%
bounds.eventgroup(2).lower = [zeros(1,7)];
bounds.eventgroup(2).upper = [zeros(1,7)];
%----------------------------------------------------------------------------%
% Event Group 3: Terminal Time of Phase 2 = 0.25*(Terminal Time of Phase 3)  %
%----------------------------------------------------------------------------%
bounds.eventgroup(3).lower = 0;
bounds.eventgroup(3).upper = 0;
%----------------------------------------------------------------------------%
% Event Group 4: Initial Value of Sixth Component of State = Parameter       %
%----------------------------------------------------------------------------%
bounds.eventgroup(4).lower = 0;
bounds.eventgroup(4).upper = 0;
%----------------------------------------------------------------------------%
% Event Group 5: Final Value of Fourth Component of State >=1                %
%----------------------------------------------------------------------------%
bounds.eventgroup(5).lower = y4fmin;
bounds.eventgroup(5).upper = y4fmax;

%-------------------------------------%
%          Setup for GPOPS-II         %
%-------------------------------------%
setup.name = 'Kinetic-Batch-Reactor';
setup.functions.continuous = @kineticBatchReactorContinuous;
setup.functions.endpoint   = @kineticBatchReactorEndpoint;

%-------------------------------------%
%      Initial Mesh for Problem       %
%-------------------------------------%
setup.mesh.phase = meshphase;

%-------------------------------------%
%      Mesh refinement options        %
%-------------------------------------%
setup.displaylevel                   = 2;
setup.mesh.method                    = 'hp-PattersonRao';
setup.mesh.tolerance                 = 1e-7;
setup.mesh.maxiterations             = 20;
setup.mesh.colpointsmin              = 3;
setup.mesh.colpointsmax              = 6;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'mumps';
setup.nlp.ipoptoptions.tolerance     = 1e-7;
setup.derivatives.supplier           = 'adigator';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';

%----------------------------------%
%   Solve Problem Using GPOPS-II   %
%----------------------------------%
output = gpops2(setup);
