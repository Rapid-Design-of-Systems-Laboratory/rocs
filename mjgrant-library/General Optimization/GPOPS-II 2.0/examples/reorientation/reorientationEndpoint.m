function output = reorientationEndpoint(input)

% Inputs
% input.phase(phasenumber).initialstate -- row
% input.phase(phasenumber).finalstate -- row
% input.phase(phasenumber).initialtime -- scalar
% input.phase(phasenumber).finaltime -- scalar
% input.phase(phasenumber).integral -- row
%
% input.parameter -- row


% Output
% output.objective -- scalar
% output.eventgroup(eventnumber).event -- row

tf = input.phase.finaltime(1);

% cost
output.objective = tf;

