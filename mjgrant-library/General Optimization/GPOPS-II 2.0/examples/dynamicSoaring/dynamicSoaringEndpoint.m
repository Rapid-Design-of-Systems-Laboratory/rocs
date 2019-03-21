function output = dynamicSoaringEndpoint(input)

t0    = input.phase(1).initialtime;
tf    = input.phase(1).finaltime;
x0    = input.phase(1).initialstate;
xf    = input.phase(1).finalstate;
beta  = input.parameter;
output.eventgroup(1).event = [xf(4)-x0(4), xf(5)-x0(5), xf(6)-x0(6)];
output.objective = 7*beta;
