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

function output = goddardRocketEndpointADiGatorGrd(input)
global ADiGator_goddardRocketEndpointADiGatorGrd
if isempty(ADiGator_goddardRocketEndpointADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_goddardRocketEndpointADiGatorGrd.goddardRocketEndpointADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
auxdata = input.auxdata;
%User Line: auxdata = input.auxdata;
%User Line: % Variables at Start and Terminus of Phase 1
cada1f1dv = input.phase(1).initialtime.dv;
cada1f1 = input.phase(1).initialtime.f;
t0{1}.dv = cada1f1dv;
t0{1}.f = cada1f1;
%User Line: t0{1} = input.phase(1).initialtime;
cada1f1dv = input.phase(1).finaltime.dv;
cada1f1 = input.phase(1).finaltime.f;
tf{1}.dv = cada1f1dv;
tf{1}.f = cada1f1;
%User Line: tf{1} = input.phase(1).finaltime;
cada1f1dv = input.phase(1).initialstate.dv;
cada1f1 = input.phase(1).initialstate.f;
x0{1}.dv = cada1f1dv;
x0{1}.f = cada1f1;
%User Line: x0{1} = input.phase(1).initialstate;
cada1f1dv = input.phase(1).finalstate.dv;
cada1f1 = input.phase(1).finalstate.f;
xf{1}.dv = cada1f1dv;
xf{1}.f = cada1f1;
%User Line: xf{1} = input.phase(1).finalstate;
%User Line: % Variables at Start and Terminus of Phase 2
cada1f1dv = input.phase(2).initialtime.dv;
cada1f1 = input.phase(2).initialtime.f;
t0{2}.dv = cada1f1dv;
t0{2}.f = cada1f1;
%User Line: t0{2} = input.phase(2).initialtime;
cada1f1dv = input.phase(2).finaltime.dv;
cada1f1 = input.phase(2).finaltime.f;
tf{2}.dv = cada1f1dv;
tf{2}.f = cada1f1;
%User Line: tf{2} = input.phase(2).finaltime;
cada1f1dv = input.phase(2).initialstate.dv;
cada1f1 = input.phase(2).initialstate.f;
x0{2}.dv = cada1f1dv;
x0{2}.f = cada1f1;
%User Line: x0{2} = input.phase(2).initialstate;
cada1f1dv = input.phase(2).finalstate.dv;
cada1f1 = input.phase(2).finalstate.f;
xf{2}.dv = cada1f1dv;
xf{2}.f = cada1f1;
%User Line: xf{2} = input.phase(2).finalstate;
%User Line: % Variables at Start and Terminus of Phase 3
cada1f1dv = input.phase(3).initialtime.dv;
cada1f1 = input.phase(3).initialtime.f;
t0{3}.dv = cada1f1dv;
t0{3}.f = cada1f1;
%User Line: t0{3} = input.phase(3).initialtime;
cada1f1dv = input.phase(3).finaltime.dv;
cada1f1 = input.phase(3).finaltime.f;
tf{3}.dv = cada1f1dv;
tf{3}.f = cada1f1;
%User Line: tf{3} = input.phase(3).finaltime;
cada1f1dv = input.phase(3).initialstate.dv;
cada1f1 = input.phase(3).initialstate.f;
x0{3}.dv = cada1f1dv;
x0{3}.f = cada1f1;
%User Line: x0{3} = input.phase(3).initialstate;
cada1f1dv = input.phase(3).finalstate.dv;
cada1f1 = input.phase(3).finalstate.f;
xf{3}.dv = cada1f1dv;
xf{3}.f = cada1f1;
%User Line: xf{3} = input.phase(3).finalstate;
%User Line: % Event Group 1:  Linkage Constraints Between Phases 1 and 2
cada1f1dv = x0{2}.dv(Gator1Data.Index2);
cada1f1 = x0{2}.f(Gator1Data.Index1);
cada1f2dv = xf{1}.dv(Gator1Data.Index4);
cada1f2 = xf{1}.f(Gator1Data.Index3);
cada1td1 = zeros(6,1);
cada1td1(Gator1Data.Index5) = cada1f1dv;
cada1td1(Gator1Data.Index6) = cada1td1(Gator1Data.Index6) + -cada1f2dv;
cada1f3dv = cada1td1;
cada1f3 = cada1f1 - cada1f2;
cada1f4dv = t0{2}.dv;
cada1f4 = t0{2}.f;
cada1f5dv = tf{1}.dv;
cada1f5 = tf{1}.f;
cada1td1 = zeros(2,1);
cada1td1(2) = cada1f4dv;
cada1td1(1) = cada1td1(1) + -cada1f5dv;
cada1f6dv = cada1td1;
cada1f6 = cada1f4 - cada1f5;
cada1td1 = zeros(8,1);
cada1td1(Gator1Data.Index7) = cada1f3dv;
cada1td1(Gator1Data.Index8) = cada1f6dv;
cada1f7dv = cada1td1;
cada1f7 = [cada1f3 cada1f6];
output.eventgroup(1).event.dv = cada1f7dv;
output.eventgroup(1).event.f = cada1f7;
%User Line: output.eventgroup(1).event = [x0{2}(1:3)-xf{1}(1:3),t0{2}-tf{1}];
%User Line: % Event Group 2:  Event that Defines Terminus of Singular Arc
h.dv = x0{2}.dv(1);
h.f = x0{2}.f(1);
%User Line: h = x0{2}(1);
v.dv = x0{2}.dv(2);
v.f = x0{2}.f(2);
%User Line: v = x0{2}(2);
m.dv = x0{2}.dv(3);
m.f = x0{2}.f(3);
%User Line: m = x0{2}(3);
cada1f1dv = 2.*v.f.^(2-1).*v.dv;
cada1f1 = v.f^2;
cada1f2dv = auxdata.dragk.*cada1f1dv;
cada1f2 = auxdata.dragk*cada1f1;
cada1f3dv = -h.dv;
cada1f3 = uminus(h.f);
cada1f4dv = cada1f3dv./auxdata.H;
cada1f4 = cada1f3/auxdata.H;
cada1f5dv = exp(cada1f4).*cada1f4dv;
cada1f5 = exp(cada1f4);
cada1td1 = zeros(2,1);
cada1td1(2) = cada1f5.*cada1f2dv;
cada1td1(1) = cada1td1(1) + cada1f2.*cada1f5dv;
D.dv = cada1td1;
D.f = cada1f2*cada1f5;
%User Line: D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
cada1f1dv = auxdata.g0.*m.dv;
cada1f1 = m.f*auxdata.g0;
cada1f2dv = v.dv./auxdata.c;
cada1f2 = v.f/auxdata.c;
cada1f3dv = cada1f2dv;
cada1f3 = 1 + cada1f2;
cada1td1 = zeros(2,1);
cada1td1(2) = D.f.*cada1f3dv;
cada1td1 = cada1td1 + cada1f3.*D.dv;
cada1f4dv = cada1td1;
cada1f4 = cada1f3*D.f;
cada1td1 = zeros(3,1);
cada1td1(3) = cada1f1dv;
cada1td1(Gator1Data.Index9) = cada1td1(Gator1Data.Index9) + -cada1f4dv;
e0.dv = cada1td1;
e0.f = cada1f1 - cada1f4;
%User Line: e0 = m*auxdata.g0-(1+v/auxdata.c).*D;
output.eventgroup(2).event.dv = e0.dv;
output.eventgroup(2).event.f = e0.f;
%User Line: output.eventgroup(2).event = e0;
%User Line: % Event Group 3:  Linkage Constraints Between Phases 2 and 3
cada1f1dv = x0{3}.dv(Gator1Data.Index11);
cada1f1 = x0{3}.f(Gator1Data.Index10);
cada1f2dv = xf{2}.dv(Gator1Data.Index13);
cada1f2 = xf{2}.f(Gator1Data.Index12);
cada1td1 = zeros(6,1);
cada1td1(Gator1Data.Index14) = cada1f1dv;
cada1td1(Gator1Data.Index15) = cada1td1(Gator1Data.Index15) + -cada1f2dv;
cada1f3dv = cada1td1;
cada1f3 = cada1f1 - cada1f2;
cada1f4dv = t0{3}.dv;
cada1f4 = t0{3}.f;
cada1f5dv = tf{2}.dv;
cada1f5 = tf{2}.f;
cada1td1 = zeros(2,1);
cada1td1(2) = cada1f4dv;
cada1td1(1) = cada1td1(1) + -cada1f5dv;
cada1f6dv = cada1td1;
cada1f6 = cada1f4 - cada1f5;
cada1td1 = zeros(8,1);
cada1td1(Gator1Data.Index16) = cada1f3dv;
cada1td1(Gator1Data.Index17) = cada1f6dv;
cada1f7dv = cada1td1;
cada1f7 = [cada1f3 cada1f6];
output.eventgroup(3).event.dv = cada1f7dv;
output.eventgroup(3).event.f = cada1f7;
%User Line: output.eventgroup(3).event = [x0{3}(1:3)-xf{2}(1:3),t0{3}-tf{2}];
%User Line: % Event Group 4:  Final Time of Each Phase Larger Than Initial Time of Phase
cada1f1dv = tf{1}.dv;
cada1f1 = tf{1}.f;
cada1f2dv = t0{1}.dv;
cada1f2 = t0{1}.f;
cada1td1 = zeros(2,1);
cada1td1(2) = cada1f1dv;
cada1td1(1) = cada1td1(1) + -cada1f2dv;
cada1f3dv = cada1td1;
cada1f3 = cada1f1 - cada1f2;
cada1f4dv = tf{2}.dv;
cada1f4 = tf{2}.f;
cada1f5dv = t0{2}.dv;
cada1f5 = t0{2}.f;
cada1td1 = zeros(2,1);
cada1td1(2) = cada1f4dv;
cada1td1(1) = cada1td1(1) + -cada1f5dv;
cada1f6dv = cada1td1;
cada1f6 = cada1f4 - cada1f5;
cada1f7dv = tf{3}.dv;
cada1f7 = tf{3}.f;
cada1f8dv = t0{3}.dv;
cada1f8 = t0{3}.f;
cada1td1 = zeros(2,1);
cada1td1(2) = cada1f7dv;
cada1td1(1) = cada1td1(1) + -cada1f8dv;
cada1f9dv = cada1td1;
cada1f9 = cada1f7 - cada1f8;
cada1td1 = zeros(6,1);
cada1td1(Gator1Data.Index18) = cada1f3dv;
cada1td1(Gator1Data.Index19) = cada1f6dv;
cada1td1(Gator1Data.Index20) = cada1f9dv;
cada1f10dv = cada1td1;
cada1f10 = [cada1f3 cada1f6 cada1f9];
output.eventgroup(4).event.dv = cada1f10dv;
output.eventgroup(4).event.f = cada1f10;
%User Line: output.eventgroup(4).event = [tf{1}-t0{1},tf{2}-t0{2},tf{3}-t0{3}];
%User Line: % Objective Function:  Maximize Final Altitude
cada1f1dv = xf{3}.dv(1);
cada1f1 = xf{3}.f(1);
output.objective.dv = -cada1f1dv;
output.objective.f = uminus(cada1f1);
%User Line: output.objective = -xf{3}(1);
%User Line: %---------------------------------------------------------------------------------%
%User Line: %--------------------- END Function goddardRocketEndpoint.m ----------------------%
%User Line: %---------------------------------------------------------------------------------%
output.eventgroup(1).event.dv_size = [4,24];
output.eventgroup(1).event.dv_location = Gator1Data.Index21;
output.eventgroup(2).event.dv_size = 24;
output.eventgroup(2).event.dv_location = Gator1Data.Index22;
output.eventgroup(3).event.dv_size = [4,24];
output.eventgroup(3).event.dv_location = Gator1Data.Index23;
output.eventgroup(4).event.dv_size = [3,24];
output.eventgroup(4).event.dv_location = Gator1Data.Index24;
output.objective.dv_size = 24;
output.objective.dv_location = Gator1Data.Index25;
end


function ADiGator_LoadData()
global ADiGator_goddardRocketEndpointADiGatorGrd
ADiGator_goddardRocketEndpointADiGatorGrd = load('goddardRocketEndpointADiGatorGrd.mat');
return
end