function output = orbitRaisingEndpoint(input)

% Inputs
% input.phase(phasenumber).initialstate -- row
% input.phase(phasenumber).finalstate -- row
% input.phase(phasenumber).initialtime -- scalar
% input.phase(phasenumber).finaltime -- scalar
% input.phase(phasenumber).integral -- row
%
% input.parameter -- row

% input.auxdata = auxiliary information
mu = input.auxdata.mu;

% Output
% output.objective -- scalar
% output.eventgroup(eventnumber).event -- row

rf = input.phase.finalstate(1);
vthetaf = input.phase.finalstate(4);

% cost
output.objective = -rf;

% event
output.eventgroup.event = vthetaf-sqrt(mu./rf);