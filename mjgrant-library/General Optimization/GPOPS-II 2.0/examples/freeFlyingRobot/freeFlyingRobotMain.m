%------------------- Free Flying Robot Problem -------------------%
% This example is taken verbatim from the following reference:    %
% Sakawa, Y, "Trajectory Planning of a Free-Flying Robot by Using %
% the Optimal Control," Optimal Control Applications and Methods, %
% Vol. 20, 1999, pp. 235-248.                                     %
%-----------------------------------------------------------------%

% ------------------------------------------------------------- %
%               Set up auxiliary data for problem               %
% ------------------------------------------------------------- %
auxdata.alpha = 0.2;
auxdata.beta  = 0.2;
auxdata.gamma = 1;

% --------------------------------------------------------------%
%           Set up bounds on state, control, and time           %
% --------------------------------------------------------------%
t0 = 0;
tf = 12;
x0       = -10;  xf       = 0;
y0       = -10;  yf       = 0;
theta0   = pi/2; thetaf   = 0;
vx0      = 0;    vxf      = 0;
vy0      = 0;    vyf      = 0;
omega0   = 0;    omegaf   = 0;

xmin     = -10; xmax     = 10;
ymin     = -10; ymax     = 10;
thetamin = -pi; thetamax = pi;
vxmin    = -2;  vxmax    = 2;
vymin    = -2;  vymax    = 2;
omegamin = -1;  omegamax = 1;

% The control is six-dimensional in this formulation.  
% The "real" controls are u1-u2 and u3-u4.  In order to
% use this formulation it is necessary to include the 
% control constraints ui>=0 (i=1,...4) along with the 
% path constraints u1+u2<=1 and u3+u4<=1.
u1Min = 0; u1Max = 1000;
u2Min = 0; u2Max = 1000;
u3Min = 0; u3Max = 1000;
u4Min = 0; u4Max = 1000;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower   = tf;
bounds.phase.finaltime.upper   = tf;

bounds.phase.initialstate.lower = [x0,y0,theta0,vx0,vy0,omega0];
bounds.phase.initialstate.upper = [x0,y0,theta0,vx0,vy0,omega0];
bounds.phase.state.lower        = [xmin,ymin,thetamin,vxmin,vymin,omegamin];
bounds.phase.state.upper        = [xmax,ymax,thetamax,vxmax,vymax,omegamax];
bounds.phase.finalstate.lower   = [xf,yf,thetaf,vxf,vyf,omegaf];
bounds.phase.finalstate.upper   = [xf,yf,thetaf,vxf,vyf,omegaf];
bounds.phase.control.lower = [u1Min,u2Min,u3Min,u4Min];
bounds.phase.control.upper = [u1Max,u2Max,u3Max,u4Max];
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 100;
bounds.phase.path.lower = [-1000, -1000];
bounds.phase.path.upper = [1, 1];

tGuess     = [t0; tf];
xGuess     = [x0; xf];
yGuess     = [y0; yf];
thetaGuess = [theta0; thetaf];
vxGuess    = [vx0; vxf];
vyGuess    = [vy0; vyf];
omegaGuess = [omega0; omegaf];
u1Guess    = [0; 0];
u2Guess    = [0; 0];
u3Guess    = [0; 0];
u4Guess    = [0; 0];
guess.phase.time = tGuess;
guess.phase.state = [xGuess,yGuess,thetaGuess,vxGuess,vyGuess,omegaGuess];
guess.phase.control = [u1Guess,u2Guess,u3Guess,u4Guess];
guess.phase.integral = 0;

mesh.method = 'hp-LiuRao-Legendre';
mesh.tolerance = 1e-6; 
mesh.maxiterations = 25;
mesh.colpointsmin = 4;
mesh.colpointsmax = 10;
mesh.phase.colpoints = 4*ones(1,10);
mesh.phase.fraction = 0.1*ones(1,10);

% ----------------------------------------------------------------------- %
% Set up solver
% ----------------------------------------------------------------------- %
setup.name = 'Free-Flying-Robot';
setup.functions.continuous = @freeFlyingRobotContinuous;
setup.functions.endpoint = @freeFlyingRobotEndpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.mesh = mesh;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.maxiterations = 3000;
setup.nlp.ipoptoptions.tolerance = 1e-7;
setup.displaylevel = 2;
setup.scales.method = 'automatic-bounds';

% ----------------------------------------------------------------------- %
% Solve problem and extract solution
% ----------------------------------------------------------------------- %
output = gpops2(setup);
