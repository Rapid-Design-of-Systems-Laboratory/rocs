%-------------------------------------------------%
% Begin Function:  minimumTimeToClimbContinuous.m %
%-------------------------------------------------%
function phaseout = brysonMinimumTimeToClimbContinuous(input)

CONSTANTS         = input.auxdata;
us1976            = CONSTANTS.us1976;
Ttab              = CONSTANTS.T;
mu                = CONSTANTS.mu;
S                 = CONSTANTS.S;
g0                = CONSTANTS.g0;
Isp               = CONSTANTS.Isp;
Re                = CONSTANTS.Re;

x                 = input.phase(1).state;
u                 = input.phase(1).control;
h                 = x(:,1);
v                 = x(:,2);
fpa               = x(:,3);
mass              = x(:,4);
alpha             = u(:,1);
r                 = h+Re;
rho               = interp1(us1976(:,1),us1976(:,2),h,'spline');
sos               = interp1(us1976(:,1),us1976(:,3),h,'spline');
Mach              = v./sos;
[CD0,Clalpha,eta] = brysonMinimumTimeToClimbCompute(Mach,CONSTANTS);
Thrust            = interp2(CONSTANTS.aa,CONSTANTS.mm,Ttab,h,Mach,'spline');
CD                = CD0 + eta.*Clalpha.*alpha.^2;
CL                = Clalpha.*alpha;
q                 = 0.5.*rho.*v.*v;
D                 = q.*S.*CD;
L                 = q.*S.*CL;
hdot              = v.*sin(fpa);
vdot              = (Thrust.*cos(alpha)-D)./mass - mu.*sin(fpa)./r.^2;
fpadot            = (Thrust.*sin(alpha)+L)./(mass.*v)+cos(fpa).*(v./r-mu./(v.*r.^2));
mdot              = -Thrust./(g0.*Isp);

phaseout.dynamics = [hdot, vdot, fpadot, mdot];

%-------------------------------------------------%
% End Function:  minimumTimeToClimbContinuous.m   %
%-------------------------------------------------%
