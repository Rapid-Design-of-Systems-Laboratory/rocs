function phaseout = freeFlyingRobotContinuous(input);

auxdata = input.auxdata;
alpha = input.auxdata.alpha;
beta  = input.auxdata.beta;

x     = input.phase.state(:,1);
y     = input.phase.state(:,2);
theta = input.phase.state(:,3);
vx    = input.phase.state(:,4);
vy    = input.phase.state(:,5);
omega = input.phase.state(:,6);
u1    = input.phase.control(:,1);
u2    = input.phase.control(:,2);
u3    = input.phase.control(:,3);
u4    = input.phase.control(:,4);
T1    = u1-u2;
T2    = u3-u4;

xdot     = vx;
ydot     = vy;
thetadot = omega;
vxdot    = (T1+T2).*cos(theta);
vydot    = (T1+T2).*sin(theta);
omegadot = alpha*T1-beta*T2;

phaseout.dynamics  = [xdot, ydot, thetadot, vxdot, vydot, omegadot];
phaseout.path      = [u1+u2, u3+u4];
phaseout.integrand = u1+u2+u3+u4;