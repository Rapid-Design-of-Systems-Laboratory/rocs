function phaseout = doubleIntegratorContinuous(input)

t  = input.phase.time;
x1 = input.phase.state(:,1);
x2 = input.phase.state(:,2);
u  = input.phase.control(:,1);

x1dot = x2;
x2dot = u;

phaseout.dynamics = [x1dot, x2dot];
