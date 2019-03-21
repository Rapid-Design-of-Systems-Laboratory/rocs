%--------------------------------------------------------------------------------%
%--------------------- BEGIN Function HMME_phases_Endpoint.m -------------------%
%--------------------------------------------------------------------------------%
function output = Thrust_real_updated_Endpoint(input)

vf = input.phase.finalstate(3);
massf = input.phase.finalstate(5);
thettaf = input.phase.finalstate(2);
tf = input.phase.finaltime;
% cost
output.objective = tf;
%output.objective = -vf;
%output.objective = -massf;
%output.objective = -thettaf;

%---------------------------------------------------------------------------------%
%--------------------- END Function HMME_phases_Endpoint.m ----------------------%
%---------------------------------------------------------------------------------%
