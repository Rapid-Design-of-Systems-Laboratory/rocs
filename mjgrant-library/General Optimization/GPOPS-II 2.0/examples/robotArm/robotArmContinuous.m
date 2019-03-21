function phaseout = moonlanderContinuous(input)

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

x = input.phase.state;
u = input.phase.control;

Iphi   = ((input.auxdata.L-x(:,1)).^3+x(:,1).^3)/3;
Itheta = Iphi.*sin(x(:,5)).^2;
x1dot = x(:,2);
x2dot = u(:,1)./input.auxdata.L;
x3dot = x(:,4);
x4dot = u(:,2)./Itheta;
x5dot = x(:,6);
x6dot = u(:,3)./Iphi;

phaseout.dynamics = [x1dot, x2dot, x3dot, x4dot, x5dot, x6dot];

