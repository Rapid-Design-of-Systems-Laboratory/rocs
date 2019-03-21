%-------------------------------------------------------------------------%
%--------------------- BEGIN Function launchEndpoint.m -------------------%
%-------------------------------------------------------------------------%
function output = launchEndpoint(input)

% Variables at Start and Terminus of Phase 1
t01 = input.phase(1).initialtime;
tf1 = input.phase(1).finaltime;
x01 = input.phase(1).initialstate;
xf1 = input.phase(1).finalstate;
% Variables at Start and Terminus of Phase 2
t02 = input.phase(2).initialtime;
tf2 = input.phase(2).finaltime;
x02 = input.phase(2).initialstate;
xf2 = input.phase(2).finalstate;
% Variables at Start and Terminus of Phase 3
t03 = input.phase(3).initialtime;
tf3 = input.phase(3).finaltime;
x03 = input.phase(3).initialstate;
xf3 = input.phase(3).finalstate;
% Variables at Start and Terminus of Phase 2
t04 = input.phase(4).initialtime;
tf4 = input.phase(4).finaltime;
x04 = input.phase(4).initialstate;
xf4 = input.phase(4).finalstate;

% Event Group 1:  Linkage Constraints Between Phases 1 and 2
output.eventgroup(1).event = [x02(1:7)-xf1(1:7), t02-tf1];
% Event Group 2:  Linkage Constraints Between Phases 2 and 3
output.eventgroup(2).event = [x03(1:7)-xf2(1:7), t03-tf2];
% Event Group 3:  Linkage Constraints Between Phases 3 and 4
output.eventgroup(3).event = [x04(1:7)-xf3(1:7), t04-tf3];
% Event Group 4:  Constraints on Terminal Orbit
orbitalElements = launchrv2oe(xf4(1:3).',xf4(4:6).',input.auxdata.mu);
output.eventgroup(4).event = orbitalElements(1:5).';
output.objective = -xf4(7);

%-------------------------------------------------------------------------%
%--------------------- END Function launchEndpoint.m ---------------------%
%-------------------------------------------------------------------------%
