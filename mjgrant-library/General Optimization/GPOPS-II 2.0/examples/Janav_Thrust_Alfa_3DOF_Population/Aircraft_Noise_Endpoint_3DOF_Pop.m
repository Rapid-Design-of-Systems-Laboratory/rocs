%--------------------------------------------------------------------------------%
%--------------------- BEGIN Function HMME_phases_Endpoint.m --------------------%
%--------------------------------------------------------------------------------%

function output = Aircraft_Noise_Endpoint_3DOF_Pop(input)
output.objective = input.phase.integral;

% tf = input.phase.finaltime;
% downf = input.phase.finalstate(2);
% output.objective = -downf; %tf;
%---------------------------------------------------------------------------------%
%--------------------- END Function HMME_phases_Endpoint.m ----------------------%
%---------------------------------------------------------------------------------%
