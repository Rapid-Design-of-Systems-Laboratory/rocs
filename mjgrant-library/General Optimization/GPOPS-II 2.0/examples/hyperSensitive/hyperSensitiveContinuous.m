%-----------------------------------%
% BEGIN: hyperSensitiveContinuous.m %
%-----------------------------------%
function phaseout = hyperSensitiveContinuous(input)


t = input.phase.time;
x = input.phase.state;
u = input.phase.control;

% xdot = -x.^3+u;
xdot = -x+u;
phaseout.dynamics = xdot;
phaseout.integrand = 0.5*(x.^2+u.^2);

%---------------------------------%
% END: hyperSensitiveContinuous.m %
%---------------------------------%
