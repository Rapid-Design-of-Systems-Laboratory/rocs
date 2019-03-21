%-------------------------------------%
% BEGIN:  kineticBatchReactorEndpoint %
%-------------------------------------%
function output = kineticBatchReactorEndpoint(input)

gamma1 = 1;
gamma2 = 100;

cmin2ms = input.auxdata.cmin2ms;
chr2min = input.auxdata.chr2min;
cmin2hr = input.auxdata.cmin2hr;
chr2sec = input.auxdata.chr2sec;
csec2hr = input.auxdata.csec2hr;
cmin2sec = input.auxdata.cmin2sec;
csec2min = input.auxdata.csec2min;

% Variables at Start and Terminus of Phase 1
t0{1} = input.phase(1).initialtime;
tf{1} = input.phase(1).finaltime;
y0{1} = input.phase(1).initialstate;
yf{1} = input.phase(1).finalstate;
p     = input.parameter(1);

% Variables at Start and Terminus of Phase 2
t0{2} = input.phase(2).initialtime;
tf{2} = input.phase(2).finaltime;
y0{2} = input.phase(2).initialstate;
yf{2} = input.phase(2).finalstate;

% Variables at Start and Terminus of Phase 3
t0{3} = input.phase(3).initialtime;
tf{3} = input.phase(3).finaltime;
y0{3} = input.phase(3).initialstate;
yf{3} = input.phase(3).finalstate;

%-----------------------------------------------------------------------------%
% Event Group 1:  Linkage Constraint on State and Time Between Phases 1 and 2 %
%-----------------------------------------------------------------------------%
output.eventgroup(1).event = [y0{2}-yf{1}, tf{1}-t0{2}];

%-----------------------------------------------------------------------------%
% Event Group 2:  Linkage Constraint on State and Time Between Phases 2 and 3 %
%-----------------------------------------------------------------------------%
output.eventgroup(2).event = [y0{3}-yf{2}, tf{2}-t0{3}];

%-----------------------------------------------------------------------------%
% Event Group 3:  Final Time of Phase 2 = 0.25 Final Time Phase 3             %
%-----------------------------------------------------------------------------%
output.eventgroup(3).event = [tf{2}-0.25*tf{3}];

%-----------------------------------------------------------------------------%
% Event Group 4:  Initial Value of y6 = Parameter                             %
%-----------------------------------------------------------------------------%
output.eventgroup(4).event = [y0{1}(6)-p];

%-----------------------------------------------------------------------------%
% Event Group 5:  Final Value y4 >= 1                                         %
%-----------------------------------------------------------------------------%
output.eventgroup(5).event = yf{3}(4);

% output.objective = gamma1*tf{3}/chr2sec+gamma2*p;
% output.objective = gamma1*tf{3}*cmin2hr+gamma2*p;
output.objective = gamma1*tf{3}+gamma2*p;

%-------------------------------------%
% BEGIN:  kineticBatchReactorEndpoint %
%-------------------------------------%
