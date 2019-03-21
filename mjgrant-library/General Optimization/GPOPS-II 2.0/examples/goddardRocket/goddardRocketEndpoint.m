%--------------------------------------------------------------------------------%
%--------------------- BEGIN Function goddardRocketEndpoint.m -------------------%
%--------------------------------------------------------------------------------%
function output = goddardRocketEndpoint(input)

auxdata = input.auxdata;

% Variables at Start and Terminus of Phase 1
t0{1} = input.phase(1).initialtime;
tf{1} = input.phase(1).finaltime;
x0{1} = input.phase(1).initialstate;
xf{1} = input.phase(1).finalstate;
% Variables at Start and Terminus of Phase 2
t0{2} = input.phase(2).initialtime;
tf{2} = input.phase(2).finaltime;
x0{2} = input.phase(2).initialstate;
xf{2} = input.phase(2).finalstate;
% Variables at Start and Terminus of Phase 3
t0{3} = input.phase(3).initialtime;
tf{3} = input.phase(3).finaltime;
x0{3} = input.phase(3).initialstate;
xf{3} = input.phase(3).finalstate;

% Event Group 1:  Linkage Constraints Between Phases 1 and 2
output.eventgroup(1).event = [x0{2}(1:3)-xf{1}(1:3),t0{2}-tf{1}];
% Event Group 2:  Event that Defines Terminus of Singular Arc
h = x0{2}(1);
v = x0{2}(2);
m = x0{2}(3);
D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
e0 = m*auxdata.g0-(1+v/auxdata.c).*D;
output.eventgroup(2).event = e0;
% Event Group 3:  Linkage Constraints Between Phases 2 and 3
output.eventgroup(3).event = [x0{3}(1:3)-xf{2}(1:3),t0{3}-tf{2}];
% Event Group 4:  Final Time of Each Phase Larger Than Initial Time of Phase
output.eventgroup(4).event = [tf{1}-t0{1},tf{2}-t0{2},tf{3}-t0{3}];

% Objective Function:  Maximize Final Altitude
output.objective = -xf{3}(1);

%---------------------------------------------------------------------------------%
%--------------------- END Function goddardRocketEndpoint.m ----------------------%
%---------------------------------------------------------------------------------%
