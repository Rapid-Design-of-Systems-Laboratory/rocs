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

g = input.auxdata.g;

t = input.phase.time;
h = input.phase.state(:,1);
v = input.phase.state(:,2);
u = input.phase.control(:,1);

dh = v;
dv = -g + u;

phaseout.dynamics = [dh, dv];
phaseout.integrand = u;
