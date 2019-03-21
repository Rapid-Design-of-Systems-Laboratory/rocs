%----------------- Maximum Range of a Hang Glider -------------------%
% This example is taken verbatim from the following reference:       %
% Bulirsh, R., Nerz, E., Pesch, H. J., and von Stryk, O., Combining  %
% Direct and Indirect Methods in Optimal Control: Range Maximization %
% of a Hang Glider, in Optimal Control, R. Bulirsch, A. Miele, J.    %
% Stoer, and K. H. Well, eds. Vol. 111, of International Series of   %
% Numerical Mathematics, Birkhauser Verlag, Basel, 1993. pp. 273-288.%
%--------------------------------------------------------------------%
clear all
clc

%-------------------------------------------------------------------%
%---------- Initialize all of the data for the problem -------------%
%-------------------------------------------------------------------%
auxdata.rho = 1.13;     % Sea Level Density (kg/m^3)
auxdata.CD0 = 0.034;    % Zero-Lift Drag Coefficient
auxdata.k   = 0.069662; % Drag Polar Parameter 
auxdata.g   = 9.80665;  % Sea Level Gravity (m/sec^2)
auxdata.m   = 100;      % mass (kg)
auxdata.S   = 14;       % Reference Area (m^2)
auxdata.uM  = 2.5;      
auxdata.R   = 100;

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
x0 = 0;   
y0  = 1000; yf = 900; 
vx0 = +13.227567500; vxf = 13.227567500; 
vy0 = -1.2875005200; vyf = -1.2875005200;

%-------------------------------------------------------------------%
%----------------------- Limits on Variables -----------------------%
%-------------------------------------------------------------------%
t0Min  = 0;   t0Max  = 0;
tfMin  = 0;   tfMax = 200;
xmin = -3000; xmax  = -xmin; 
ymin = 0;     ymax  = 3000;
vxmin = -15;  vxmax = 15; 
vymin = -15;  vymax = 15;
CLmin = 0;    CLmax  = 1.4;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase.initialtime.lower  = t0Min;
bounds.phase.initialtime.upper  = t0Max;
bounds.phase.finaltime.lower    = tfMin;
bounds.phase.finaltime.upper    = tfMax;
bounds.phase.initialstate.lower = [x0, y0, vx0, vy0];
bounds.phase.initialstate.upper = [x0, y0, vx0, vy0];
bounds.phase.state.lower        = [xmin, ymin, vxmin, vymin];
bounds.phase.state.upper        = [xmax, ymax, vxmax, vymax];
bounds.phase.finalstate.lower   = [xmin, yf, vxf, vyf];
bounds.phase.finalstate.upper   = [xmax, yf, vxf, vyf];
bounds.phase.control.lower      = CLmin;
bounds.phase.control.upper      = CLmax;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
N                    = 10;
mesh.method          = 'hp-DarbyRao';
mesh.tolerance       = 1e-6;
mesh.phase.fraction  = (1/N)*ones(1,N);
mesh.phase.colpoints = 3*ones(1,N);
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 16;
mesh.maxiterations   = 25;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess              = [0;   100];
xGuess              = [x0;  1250];
yGuess              = [y0;  yf];
vxGuess             = [vx0; vxf];
vyGuess             = [vy0; vyf];
CLGuess             = [1; 1];
guess.phase.time    = tGuess;
guess.phase.state   = [xGuess, yGuess, vxGuess, vyGuess];
guess.phase.control = CLGuess;
guess.phase.integral = 0;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the Information Provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Hang-Glider';
setup.functions.continuous           = @hangGliderContinuous;
setup.functions.endpoint             = @hangGliderEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-7;
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
setup.mesh                           = mesh;

%-------------------------------------------------------------------%
%------------------ Solve Problem Using GPOPS-II -------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
