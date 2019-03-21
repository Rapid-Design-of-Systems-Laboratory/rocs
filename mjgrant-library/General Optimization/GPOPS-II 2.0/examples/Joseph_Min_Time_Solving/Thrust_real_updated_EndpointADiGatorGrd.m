% This code was generated using ADiGator version 1.3
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function output = Thrust_real_updated_EndpointADiGatorGrd(input)
global ADiGator_Thrust_real_updated_EndpointADiGatorGrd
if isempty(ADiGator_Thrust_real_updated_EndpointADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_Thrust_real_updated_EndpointADiGatorGrd.Thrust_real_updated_EndpointADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
vf.dv = input.phase.finalstate.dv(3);
vf.f = input.phase.finalstate.f(3);
%User Line: vf = input.phase.finalstate(3);
massf.dv = input.phase.finalstate.dv(5);
massf.f = input.phase.finalstate.f(5);
%User Line: massf = input.phase.finalstate(5);
thettaf.dv = input.phase.finalstate.dv(2);
thettaf.f = input.phase.finalstate.f(2);
%User Line: thettaf = input.phase.finalstate(2);
tf.dv = input.phase.finaltime.dv; tf.f = input.phase.finaltime.f;
%User Line: tf = input.phase.finaltime;
%User Line: % cost
output.objective.dv = tf.dv; output.objective.f = tf.f;
%User Line: output.objective = tf;
%User Line: %output.objective = -vf;
%User Line: %output.objective = -massf;
%User Line: %output.objective = -thettaf;
%User Line: %---------------------------------------------------------------------------------%
%User Line: %--------------------- END Function HMME_phases_Endpoint.m ----------------------%
%User Line: %---------------------------------------------------------------------------------%
output.objective.dv_size = 14;
output.objective.dv_location = Gator1Data.Index1;
end


function ADiGator_LoadData()
global ADiGator_Thrust_real_updated_EndpointADiGatorGrd
ADiGator_Thrust_real_updated_EndpointADiGatorGrd = load('Thrust_real_updated_EndpointADiGatorGrd.mat');
return
end