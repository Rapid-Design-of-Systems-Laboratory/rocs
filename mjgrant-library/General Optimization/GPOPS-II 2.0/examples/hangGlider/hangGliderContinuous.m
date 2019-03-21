function phaseout = hangGliderContinuous(input)

uM    = input.auxdata.uM;
R     = input.auxdata.R;
m     = input.auxdata.m;
CD0   = input.auxdata.CD0;
k     = input.auxdata.k;
rho   = input.auxdata.rho;
g     = input.auxdata.g;
S     = input.auxdata.S;
m     = input.auxdata.m;

x     = input.phase.state(:,1);
y     = input.phase.state(:,2);
vx    = input.phase.state(:,3);
vy    = input.phase.state(:,4);
CL    = input.phase.control(:,1);

CD    = CD0+k.*CL.^2;
e     = x./R;
X     = (e-2.5).^2;
uA    = uM.*(1-X).*exp(-X);
Vy    = vy-uA;
vr    = sqrt(vx.^2+Vy.^2);
q     = 0.5*rho.*vr.^2;
D     = q.*S.*CD/m; % 0.5.*CD.*rho.*S.*vr.^2;
L     = q.*S.*CL/m; % 0.5.*CL.*rho.*S.*vr.^2;
sinN  = Vy./vr;
cosN  = vx./vr;
W     = m.*g;

xdot  = vx;
ydot  = vy;
vxdot = -L.*sinN-D.*cosN;
vydot = L.*cosN-D.*sinN-g;

phaseout.dynamics = [xdot, ydot, vxdot, vydot];
