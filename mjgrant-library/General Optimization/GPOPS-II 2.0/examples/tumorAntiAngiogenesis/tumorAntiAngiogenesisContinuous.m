function phaseout = tumorAntiAngiogenesisContinuous(input)

zeta = input.auxdata.zeta;
b    = input.auxdata.b;
mu   = input.auxdata.mu;
d    = input.auxdata.d;
G    = input.auxdata.G;

t = input.phase(1).time;
x = input.phase(1).state;
u = input.phase(1).control;

p = x(:,1);
q = x(:,2);
pdot = -zeta.*p.*log(p./q);
qdot = q.*(b-(mu+(d*(p.^(2/3)))+G.*u));

phaseout.dynamics = [pdot qdot];
phaseout.integrand = u;
