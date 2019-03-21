function output = goddardRocketContinuous(input)
%--------------------------------------%
% Begin File:  goddardRocketContinuous %
%--------------------------------------%

auxdata = input.auxdata;

t = input.phase(1).time;
x = input.phase(1).state;
T = input.phase(1).control;
h = x(:,1);
v = x(:,2);
m = x(:,3);
D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
hdot = v;
vdot = (T-D)./m-auxdata.g0*ones(size(t));
mdot = -T./auxdata.c;
output(1).dynamics = [hdot, vdot, mdot];

t = input.phase(2).time;
x = input.phase(2).state;
T = input.phase(2).control;
h = x(:,1);
v = x(:,2);
m = x(:,3);
D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
hdot = v;
vdot = (T-D)./m-auxdata.g0*ones(size(t));
mdot = -T./auxdata.c;
output(2).dynamics = [hdot, vdot, mdot];
voverc = v/auxdata.c;
xmg = m*auxdata.g0;
term1 = (auxdata.c^2).*(ones(size(t))+voverc)./(auxdata.H*auxdata.g0)-ones(size(t))-2./voverc;
term2 = xmg./(ones(size(t))+4./voverc+2./(voverc.^2));
output(2).path = T-D-xmg-term1.*term2;

t = input.phase(3).time;
x = input.phase(3).state;
T = input.phase(3).control;
h = x(:,1);
v = x(:,2);
m = x(:,3);
D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
hdot = v;
vdot = (T-D)./m-auxdata.g0*ones(size(t));
mdot = -T./auxdata.c;
output(3).dynamics = [hdot, vdot, mdot];

%------------------------------------%
% End File:  goddardRocketContinuous %
%------------------------------------%
