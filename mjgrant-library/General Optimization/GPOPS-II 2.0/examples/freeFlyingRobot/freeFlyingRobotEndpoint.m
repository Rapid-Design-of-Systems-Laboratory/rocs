function output = freeFlyingRobotEndpoint(input);

gamma = input.auxdata.gamma;
q = input.phase.integral;
output.objective = gamma*q;