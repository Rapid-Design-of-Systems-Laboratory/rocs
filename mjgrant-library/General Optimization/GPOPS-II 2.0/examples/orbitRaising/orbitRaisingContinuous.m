function phaseout = orbitRaisingContinuous(input)

% input
% input.phase(phasenumber).state
% input.phase(phasenumber).control
% input.phase(phasenumber).time
% input.phase(phasenumber).parameter
%
% input.auxdata = auxiliary information
%
% output
% phaseout(phasenumber).dynamics
% phaseout(phasenumber).path
% phaseout(phasenumber).integrand
 
T = input.auxdata.T;
m0 = input.auxdata.m0;
dm = input.auxdata.dm;
mu = input.auxdata.mu;

t      = input.phase.time;
r      = input.phase.state(:,1);
theta  = input.phase.state(:,2);
vr     = input.phase.state(:,3);
vtheta = input.phase.state(:,4);
u1 = input.phase.control(:,1);
u2 = input.phase.control(:,2);

a       = T./(m0 - dm.*t);
dr      = vr;
dtheta  = vtheta./r;
dvr     = (vtheta.^2)./r-mu./(r.^2)+a.*u1;
dvtheta = -(vr.*vtheta)./r+a.*u2;

phaseout.dynamics  = [dr, dtheta, dvr, dvtheta];
phaseout.path = u1.^2 + u2.^2;
