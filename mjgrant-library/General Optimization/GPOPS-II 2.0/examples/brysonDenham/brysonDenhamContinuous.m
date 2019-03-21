function phaseout = brysonDenhamContinuous(input)

v                  = input.phase.state(:,2);
u                  = input.phase.control;
xdot               = v;
vdot               = u;
integrand          = (u.^2)./2;
phaseout.dynamics  = [xdot,vdot];
phaseout.integrand = integrand;
phaseout.path      = input.phase.state(:,1);
