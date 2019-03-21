%--------- Minimum Time-to-Climb of a Supersonic Aircraft ----------%
% This example is taken verbatim from the following reference:      %
% Bryson, A. E., Desai, M. N. and Hoffman, W. C., "Energy-State     %
% Approximation in Performance Optimization of Supersonic           %
% Aircraft," Journal of Aircraft, Vol. 6, No. 6, November-December, %
% 1969, pp. 481-488.                                                %
%-------------------------------------------------------------------%
clear all
close all
clc

%-------------------------------------------------------------------%
%---------- Initialize all of the data for the problem -------------%
%-------------------------------------------------------------------%
load brysonMinimumTimeToClimbAeroData.mat;

%-------------------------------------------------------------------%
%----------------- U.S. 1976 Standard Atmosphere  ------------------%
%-------------------------------------------------------------------%
% Format of Data:                                                   %
%   Column 1:  Altitude (m)                                         %
%   Column 2:  Atmospheric Density (kg/m^3)                         %
%   Column 3:  Speed of Sound (m/s)                                 %
%-------------------------------------------------------------------%
us1976 = [-2000    1.478e+00     3.479e+02
              0    1.225e+00     3.403e+02
           2000    1.007e+00     3.325e+02
           4000    8.193e-01     3.246e+02
           6000    6.601e-01     3.165e+02
           8000    5.258e-01     3.081e+02
          10000    4.135e-01     2.995e+02
          12000    3.119e-01     2.951e+02
          14000    2.279e-01     2.951e+02
          16000    1.665e-01     2.951e+02
          18000    1.216e-01     2.951e+02
          20000    8.891e-02     2.951e+02
          22000    6.451e-02     2.964e+02
          24000    4.694e-02     2.977e+02
          26000    3.426e-02     2.991e+02
          28000    2.508e-02     3.004e+02
          30000    1.841e-02     3.017e+02
          32000    1.355e-02     3.030e+02
          34000    9.887e-03     3.065e+02
          36000    7.257e-03     3.101e+02
          38000    5.366e-03     3.137e+02
          40000    3.995e-03     3.172e+02
          42000    2.995e-03     3.207e+02
          44000    2.259e-03     3.241e+02
          46000    1.714e-03     3.275e+02
          48000    1.317e-03     3.298e+02
          50000    1.027e-03     3.298e+02
          52000    8.055e-04     3.288e+02
          54000    6.389e-04     3.254e+02
          56000    5.044e-04     3.220e+02
          58000    3.962e-04     3.186e+02
          60000    3.096e-04     3.151e+02
          62000    2.407e-04     3.115e+02
          64000    1.860e-04     3.080e+02
          66000    1.429e-04     3.044e+02
          68000    1.091e-04     3.007e+02
          70000    8.281e-05     2.971e+02
          72000    6.236e-05     2.934e+02
          74000    4.637e-05     2.907e+02
          76000    3.430e-05     2.880e+02
          78000    2.523e-05     2.853e+02
          80000    1.845e-05     2.825e+02
          82000    1.341e-05     2.797e+02
          84000    9.690e-06     2.769e+02
          86000    6.955e-06     2.741e+02];

%-------------------------------------------------------------------%
%--------------- Propulsion Data for Bryson Aircraft ---------------%
%-------------------------------------------------------------------%
% The thrust depends for the aircraft considered by Bryson in 1969  %
% depends upon the Mach number and the altitude.  This data is taken%
% verbatim from the 1969 Journal of Aircraft paper (see reference   %
% above) and is copied for use in this example.  The data are stored%
% in the following variables:                                       %
%    - Mtab: a vector of values of Mach number                      %
%    - alttab: a vector of altitude values                          %
%    - Ttab: is a table of aircraft thrust values                   %
% After conversion, the altitude given in meters.                   %
% After conversion, the thrust given in Newtons.                    %
%-------------------------------------------------------------------%
Mtab   = [0; 0.2; 0.4; 0.6; 0.8; 1; 1.2; 1.4; 1.6; 1.8];
alttab = 304.8*[0 5 10 15 20 25 30 40 50 70];
Ttab   = 4448.222*[24.2 24.0 20.3 17.3 14.5 12.2 10.2  5.7  3.4 0.1;
                   28.0 24.6 21.1 18.1 15.2 12.8 10.7  6.5  3.9 0.2;
                   28.3 25.2 21.9 18.7 15.9 13.4 11.2  7.3  4.4 0.4;
                   30.8 27.2 23.8 20.5 17.3 14.7 12.3  8.1  4.9 0.8;
                   34.5 30.3 26.6 23.2 19.8 16.8 14.1  9.4  5.6 1.1;
                   37.9 34.3 30.4 26.8 23.3 19.8 16.8 11.2  6.8 1.4;
                   36.1 38.0 34.9 31.3 27.3 23.6 20.1 13.4  8.3 1.7;
                   36.1 36.6 38.5 36.1 31.6 28.1 24.2 16.2 10.0 2.2;
                   36.1 35.2 42.1 38.7 35.7 32.0 28.1 19.3 11.9 2.9;
                   36.1 33.8 45.7 41.3 39.8 34.6 31.1 21.7 13.3 3.1];

%-------------------------------------------------------------------%
%-------------- Aerodynamic Data for Bryson Aircraft ---------------%
%-------------------------------------------------------------------%
%  M2 is a vector of Mach number values                             %
%  Clalphatab is a vector of coefficient of lift values             %
%  CD0tab is a vector of zero-lift coefficient of drag values       %
%   - etatab is a vector of load factors                            %
%-------------------------------------------------------------------%
M2         = [0 0.4 0.8 0.9 1.0 1.2 1.4 1.6 1.8];
Clalphatab = [3.44 3.44 3.44 3.58 4.44 3.44 3.01 2.86 2.44];
CD0tab     = [0.013 0.013 0.013 0.014 0.031 0.041 0.039 0.036 0.035];
etatab     = [0.54 0.54 0.54 0.75 0.79 0.78 0.89 0.93 0.93];


%-------------------------------------------------------------------%
%---- All Data Required by User Functions is Stored in AUXDATA -----%
%-------------------------------------------------------------------%
auxdata.CDdat     = CDdat;
auxdata.CLdat     = CLdat;
auxdata.etadat    = etadat;
auxdata.M         = Mtab;
auxdata.M2        = M2;
auxdata.alt       = alttab;
auxdata.T         = Ttab;
auxdata.Clalpha   = Clalphatab;
auxdata.CD0       = CD0tab;
auxdata.eta       = etatab;
auxdata.ppCLalpha = polyfit(auxdata.M2,auxdata.Clalpha,8);
auxdata.ppCD0     = polyfit(auxdata.M2,auxdata.CD0,8);
auxdata.ppeta     = polyfit(auxdata.M2,auxdata.eta,8);
auxdata.Re        = 6378145;
auxdata.mu        = 3.986e14;
auxdata.S         = 49.2386;
auxdata.g0        = 9.80665;
auxdata.Isp       = 1600;
auxdata.H         = 7254.24;
auxdata.rho0      = 1.225;
auxdata.us1976    = us1976;
[aa,mm]           = meshgrid(alttab,Mtab);
auxdata.aa        = aa;
auxdata.mm        = mm;

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
t0     = 0;         % Initial time (sec)
alt0   = 0;         % Initial altitude (meters)
altf   = 19994.88;  % Final altitude (meters)
speed0 = 129.314;   % Initial speed (m/s)
speedf = 295.092;   % Final speed (m/s)
fpa0   = 0;         % Initial flight path angle (rad)
fpaf   = 0;         % Final flight path angle (rad)
mass0  = 19050.864; % Initial mass (kg)

%-------------------------------------------------------------------%
%----------------------- Limits on Variables -----------------------%
%-------------------------------------------------------------------%
tfmin    = 100;        tfmax    = 800;
altmin   = 0;          altmax   = 21031.2;
speedmin = 5;          speedmax = 1000;
fpamin   = -40*pi/180; fpamax   = 40*pi/180;
massmin  = 22;         massmax  = 20410;
alphamin = -pi/4;      alphamax = pi/4;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase.initialtime.lower  = t0;
bounds.phase.initialtime.upper  = t0;
bounds.phase.finaltime.lower    = tfmin;
bounds.phase.finaltime.upper    = tfmax;
bounds.phase.initialstate.lower = [alt0, speed0, fpa0, mass0];
bounds.phase.initialstate.upper = [alt0, speed0, fpa0, mass0];
bounds.phase.state.lower        = [altmin, speedmin, fpamin, massmin];
bounds.phase.state.upper        = [altmax, speedmax, fpamax, massmax];
bounds.phase.finalstate.lower   = [altf, speedf, fpaf, massmin];
bounds.phase.finalstate.upper   = [altf, speedf, fpaf, massmax];
bounds.phase.control.lower      = alphamin;
bounds.phase.control.upper      = alphamax;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time                = [0; 1000];
guess.phase.state(:,1)          = [alt0; altf];
guess.phase.state(:,2)          = [speed0; speedf];
guess.phase.state(:,3)          = [fpa0; fpaf];
guess.phase.state(:,4)          = [mass0; mass0];
guess.phase.control             = [20; -20]*pi/180;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
mesh.tolerance    = 1e-6;
mesh.colpointsmin = 4;
mesh.colpointsmax = 10;
mesh.sigma        = 0.75;

% mesh.colpointsmax = 100;
% mesh.colpoints    = 10;
% mesh.fraction     = 1;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Bryson-Minimum-Time-to-Climb-Problem';
setup.functions.continuous           = @brysonMinimumTimeToClimbContinuous;
setup.functions.endpoint             = @brysonMinimumTimeToClimbEndpoint;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.auxdata                        = auxdata;
setup.derivatives.supplier           = 'sparseCD';
% setup.derivatives.supplier           = 'adigator';
setup.derivatives.derivativelevel    = 'second';
setup.derivatives.dependencies       = 'sparseNaN';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);

