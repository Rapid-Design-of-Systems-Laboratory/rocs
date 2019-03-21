function output = lowThrustEndpoint(input)

ff = input.phase.finalstate(2);
gf = input.phase.finalstate(3);
hf = input.phase.finalstate(4);
kf = input.phase.finalstate(5);

wf = input.phase.finalstate(7);

% Declare objective (minimize final weight)
output.objective = -wf;

% Declare event constraints
output.eventgroup.event = [ff^2+gf^2 hf^2+kf^2 ff*hf+gf*kf gf*hf-kf*ff];
