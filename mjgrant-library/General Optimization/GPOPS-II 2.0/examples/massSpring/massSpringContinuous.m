%-----------------------------------%
% BEGIN: hyperSensitiveContinuous.m %
%-----------------------------------%
function phaseout = hyperSensitiveContinuous(input)


t = input.phase.time;
x = input.phase.state;
u = input.phase.control;

xdot = u;
phaseout.dynamics = xdot;
phaseout.integrand = 0.5*(u.^2-x.^2);

%---------------------------------%
% END: hyperSensitiveContinuous.m %
%---------------------------------%
