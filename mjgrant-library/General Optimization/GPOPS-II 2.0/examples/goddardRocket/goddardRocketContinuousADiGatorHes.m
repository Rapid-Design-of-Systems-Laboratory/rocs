% This code was generated using ADiGator version 1.3
% �2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function output = goddardRocketContinuousADiGatorHes(input)
global ADiGator_goddardRocketContinuousADiGatorHes
if isempty(ADiGator_goddardRocketContinuousADiGatorHes); ADiGator_LoadData(); end
Gator1Data = ADiGator_goddardRocketContinuousADiGatorHes.goddardRocketContinuousADiGatorHes.Gator1Data;
Gator2Data = ADiGator_goddardRocketContinuousADiGatorHes.goddardRocketContinuousADiGatorHes.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: %--------------------------------------%
%User Line: % Begin File:  goddardRocketContinuous %
%User Line: %--------------------------------------%
auxdata = input.auxdata;
%User Line: auxdata = input.auxdata;
t.dV = input.phase(1).time.dV;
t.f = input.phase(1).time.f;
%User Line: t = input.phase(1).time;
x.dV = input.phase(1).state.dV;
x.f = input.phase(1).state.f;
%User Line: x = input.phase(1).state;
T.dV = input.phase(1).control.dV;
T.f = input.phase(1).control.f;
%User Line: T = input.phase(1).control;
h.dV = x.dV(:,1);
h.f = x.f(:,1);
%User Line: h = x(:,1);
v.dV = x.dV(:,2);
v.f = x.f(:,2);
%User Line: v = x(:,2);
m.dV = x.dV(:,3);
m.f = x.f(:,3);
%User Line: m = x(:,3);
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f1dVdV = v.dV.*cada2f2dV;
cada1f1dV = cada2f2.*v.dV;
cada1f1 = v.f.^2;
cada1f2dVdV = auxdata.dragk.*cada1f1dVdV;
cada1f2dV = auxdata.dragk*cada1f1dV;
cada1f2 = auxdata.dragk*cada1f1;
cada1f3dV = uminus(h.dV);
cada1f3 = uminus(h.f);
cada1f4dV = cada1f3dV/auxdata.H;
cada1f4 = cada1f3/auxdata.H;
cada2f1dV = exp(cada1f4).*cada1f4dV;
cada2f1 = exp(cada1f4);
cada1f5dVdV = cada1f4dV.*cada2f1dV;
cada1f5dV = cada2f1.*cada1f4dV;
cada1f5 = exp(cada1f4);
cada2f1 = size(cada1f2dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f5dV,1),2);
cada2td1(:,1) = cada1f2dV.*cada1f5dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f5.*cada1f2dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f5.*cada1f2dV;
cada1td1dV = cada2f1dV;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2td1 = zeros(size(cada1f2dV,1),2);
cada2td1(:,2) = cada1f5dV.*cada1f2dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f2.*cada1f5dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f2.*cada1f5dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index1) = cada2f3dV;
cada2td1(:,Gator2Data.Index2) = cada1td1dV(:,Gator2Data.Index3);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f3;
D.dVdV = cada1td1dV; D.dV = cada1td1;
D.f = cada1f2.*cada1f5;
%User Line: D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
hdot.dV = v.dV;
hdot.f = v.f;
%User Line: hdot = v;
cada2f1 = size(T.dV,1);
cada1td1 = zeros(cada2f1,3);
cada1td1(:,3) = T.dV;
cada2f1 = cada1td1(:,Gator1Data.Index1);
cada2f2dV = -D.dVdV;
cada2f2 = uminus(D.dV);
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada1td1dV = cada2f3dV;
cada1td1(:,Gator1Data.Index1) = cada2f3;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = T.f - D.f;
cada1tf1dV = m.dV(:,Gator2Data.Index4);
cada1tf1 = m.f(:,Gator1Data.Index2);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,4);
cada2tf1 = cada1tf1(:,Gator2Data.Index5);
cada2td1 = zeros(size(cada1f1dVdV,1),7);
cada2td1(:,Gator2Data.Index6) = cada1f1dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index7) = cada2td1(:,Gator2Data.Index7) + -cada1f1dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f1dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index3) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2f2dV = -cada1f1dV;
cada2f2 = uminus(cada1f1);
cada2f3dV = 2.*m.f.^(2-1).*m.dV;
cada2f3 = m.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index8);
cada2td1 = zeros(size(cada2f2dV,1),4);
cada2td1(:,Gator2Data.Index9) = cada2f2dV./cada2tf1;
cada2td1(:,3) = cada2td1(:,3) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = m.dV(:,Gator2Data.Index10);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*m.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),11);
cada2td1(:,Gator2Data.Index11) = cada2f6dV;
cada2td1(:,Gator2Data.Index12) = cada1td1dV(:,Gator2Data.Index13);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f6;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = cada1f1./m.f;
cada1f3 = size(t.f);
cada1f4 = ones(cada1f3);
cada1f5 = auxdata.g0*cada1f4;
vdot.dVdV = cada1f2dVdV; vdot.dV = cada1f2dV;
vdot.f = cada1f2 - cada1f5;
%User Line: vdot = (T-D)./m-auxdata.g0*ones(size(t));
cada1f1dV = uminus(T.dV);
cada1f1 = uminus(T.f);
mdot.dV = cada1f1dV/auxdata.c;
mdot.f = cada1f1/auxdata.c;
%User Line: mdot = -T./auxdata.c;
cada2f1 = size(hdot.f,1);
cada1td1 = zeros(cada2f1,6);
cada1td1(:,2) = hdot.dV;
cada1td1dV = vdot.dVdV;
cada1td1(:,Gator1Data.Index4) = vdot.dV;
cada1td1(:,6) = mdot.dV;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = [hdot.f vdot.f mdot.f];
output(1).dynamics.dVdV = cada1f1dVdV;
output(1).dynamics.dV = cada1f1dV;
output(1).dynamics.f = cada1f1;
%User Line: output(1).dynamics = [hdot, vdot, mdot];
t.dV = input.phase(2).time.dV;
t.f = input.phase(2).time.f;
%User Line: t = input.phase(2).time;
x.dV = input.phase(2).state.dV;
x.f = input.phase(2).state.f;
%User Line: x = input.phase(2).state;
T.dV = input.phase(2).control.dV;
T.f = input.phase(2).control.f;
%User Line: T = input.phase(2).control;
h.dV = x.dV(:,1);
h.f = x.f(:,1);
%User Line: h = x(:,1);
v.dV = x.dV(:,2);
v.f = x.f(:,2);
%User Line: v = x(:,2);
m.dV = x.dV(:,3);
m.f = x.f(:,3);
%User Line: m = x(:,3);
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f1dVdV = v.dV.*cada2f2dV;
cada1f1dV = cada2f2.*v.dV;
cada1f1 = v.f.^2;
cada1f2dVdV = auxdata.dragk.*cada1f1dVdV;
cada1f2dV = auxdata.dragk*cada1f1dV;
cada1f2 = auxdata.dragk*cada1f1;
cada1f3dV = uminus(h.dV);
cada1f3 = uminus(h.f);
cada1f4dV = cada1f3dV/auxdata.H;
cada1f4 = cada1f3/auxdata.H;
cada2f1dV = exp(cada1f4).*cada1f4dV;
cada2f1 = exp(cada1f4);
cada1f5dVdV = cada1f4dV.*cada2f1dV;
cada1f5dV = cada2f1.*cada1f4dV;
cada1f5 = exp(cada1f4);
cada2f1 = size(cada1f2dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f5dV,1),2);
cada2td1(:,1) = cada1f2dV.*cada1f5dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f5.*cada1f2dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f5.*cada1f2dV;
cada1td1dV = cada2f1dV;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2td1 = zeros(size(cada1f2dV,1),2);
cada2td1(:,2) = cada1f5dV.*cada1f2dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f2.*cada1f5dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f2.*cada1f5dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index14) = cada2f3dV;
cada2td1(:,Gator2Data.Index15) = cada1td1dV(:,Gator2Data.Index16);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f3;
D.dVdV = cada1td1dV; D.dV = cada1td1;
D.f = cada1f2.*cada1f5;
%User Line: D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
hdot.dV = v.dV;
hdot.f = v.f;
%User Line: hdot = v;
cada2f1 = size(T.dV,1);
cada1td1 = zeros(cada2f1,3);
cada1td1(:,3) = T.dV;
cada2f1 = cada1td1(:,Gator1Data.Index5);
cada2f2dV = -D.dVdV;
cada2f2 = uminus(D.dV);
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada1td1dV = cada2f3dV;
cada1td1(:,Gator1Data.Index5) = cada2f3;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = T.f - D.f;
cada1tf1dV = m.dV(:,Gator2Data.Index17);
cada1tf1 = m.f(:,Gator1Data.Index6);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,4);
cada2tf1 = cada1tf1(:,Gator2Data.Index18);
cada2td1 = zeros(size(cada1f1dVdV,1),7);
cada2td1(:,Gator2Data.Index19) = cada1f1dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index20) = cada2td1(:,Gator2Data.Index20) + -cada1f1dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f1dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index7) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2f2dV = -cada1f1dV;
cada2f2 = uminus(cada1f1);
cada2f3dV = 2.*m.f.^(2-1).*m.dV;
cada2f3 = m.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index21);
cada2td1 = zeros(size(cada2f2dV,1),4);
cada2td1(:,Gator2Data.Index22) = cada2f2dV./cada2tf1;
cada2td1(:,3) = cada2td1(:,3) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = m.dV(:,Gator2Data.Index23);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*m.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),11);
cada2td1(:,Gator2Data.Index24) = cada2f6dV;
cada2td1(:,Gator2Data.Index25) = cada1td1dV(:,Gator2Data.Index26);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f6;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = cada1f1./m.f;
cada1f3 = size(t.f);
cada1f4 = ones(cada1f3);
cada1f5 = auxdata.g0*cada1f4;
vdot.dVdV = cada1f2dVdV; vdot.dV = cada1f2dV;
vdot.f = cada1f2 - cada1f5;
%User Line: vdot = (T-D)./m-auxdata.g0*ones(size(t));
cada1f1dV = uminus(T.dV);
cada1f1 = uminus(T.f);
mdot.dV = cada1f1dV/auxdata.c;
mdot.f = cada1f1/auxdata.c;
%User Line: mdot = -T./auxdata.c;
cada2f1 = size(hdot.f,1);
cada1td1 = zeros(cada2f1,6);
cada1td1(:,2) = hdot.dV;
cada1td1dV = vdot.dVdV;
cada1td1(:,Gator1Data.Index8) = vdot.dV;
cada1td1(:,6) = mdot.dV;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = [hdot.f vdot.f mdot.f];
output(2).dynamics.dVdV = cada1f1dVdV;
output(2).dynamics.dV = cada1f1dV;
output(2).dynamics.f = cada1f1;
%User Line: output(2).dynamics = [hdot, vdot, mdot];
voverc.dV = v.dV/auxdata.c;
voverc.f = v.f/auxdata.c;
%User Line: voverc = v/auxdata.c;
xmg.dV = auxdata.g0*m.dV;
xmg.f = m.f*auxdata.g0;
%User Line: xmg = m*auxdata.g0;
cada1f1 = auxdata.c^2;
cada1f2 = size(t.f);
cada1f3 = ones(cada1f2);
cada1f4dV = voverc.dV;
cada1f4 = cada1f3 + voverc.f;
cada1f5dV = cada1f1*cada1f4dV;
cada1f5 = cada1f1*cada1f4;
cada1f6 = auxdata.H*auxdata.g0;
cada1f7dV = cada1f5dV/cada1f6;
cada1f7 = cada1f5/cada1f6;
cada1f8 = size(t.f);
cada1f9 = ones(cada1f8);
cada1f10dV = cada1f7dV;
cada1f10 = cada1f7 - cada1f9;
cada2f1dV = 2.*voverc.f.^(2-1).*voverc.dV;
cada2f1 = voverc.f.^2;
cada2f2dV = --2./cada2f1.^2.*cada2f1dV;
cada2f2 = -2./cada2f1;
cada1f11dVdV = voverc.dV.*cada2f2dV;
cada1f11dV = cada2f2.*voverc.dV;
cada1f11 = 2./voverc.f;
cada1td1 = cada1f10dV;
cada2f1dV = -cada1f11dVdV;
cada2f1 = uminus(cada1f11dV);
cada1td1dV = cada2f1dV;
cada1td1 = cada1td1 + cada2f1;
term1.dVdV = cada1td1dV; term1.dV = cada1td1;
term1.f = cada1f10 - cada1f11;
%User Line: term1 = (auxdata.c^2).*(ones(size(t))+voverc)./(auxdata.H*auxdata.g0)-ones(size(t))-2./voverc;
cada1f1 = size(t.f);
cada1f2 = ones(cada1f1);
cada2f1dV = 2.*voverc.f.^(2-1).*voverc.dV;
cada2f1 = voverc.f.^2;
cada2f2dV = --4./cada2f1.^2.*cada2f1dV;
cada2f2 = -4./cada2f1;
cada1f3dVdV = voverc.dV.*cada2f2dV;
cada1f3dV = cada2f2.*voverc.dV;
cada1f3 = 4./voverc.f;
cada1f4dVdV = cada1f3dVdV; cada1f4dV = cada1f3dV;
cada1f4 = cada1f2 + cada1f3;
cada2f1dV = 1.*voverc.f.^(1-1).*voverc.dV;
cada2f1 = voverc.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f5dVdV = voverc.dV.*cada2f2dV;
cada1f5dV = cada2f2.*voverc.dV;
cada1f5 = voverc.f.^2;
cada2f1dV = 2.*cada1f5.^(2-1).*cada1f5dV;
cada2f1 = cada1f5.^2;
cada2f2dV = --2./cada2f1.^2.*cada2f1dV;
cada2f2 = -2./cada2f1;
cada2td1 = cada1f5dV.*cada2f2dV;
cada2td1 = cada2td1 + cada2f2.*cada1f5dVdV;
cada1f6dVdV = cada2td1;
cada1f6dV = cada2f2.*cada1f5dV;
cada1f6 = 2./cada1f5;
cada1td1dV = cada1f4dVdV; cada1td1 = cada1f4dV;
cada2td1 = cada1td1dV;
cada2td1 = cada2td1 + cada1f6dVdV;
cada1td1dV = cada2td1;
cada1td1 = cada1td1 + cada1f6dV;
cada1f7dVdV = cada1td1dV; cada1f7dV = cada1td1;
cada1f7 = cada1f4 + cada1f6;
cada2f1 = size(xmg.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dV = -xmg.dV./cada1f7.^2.*cada1f7dV;
cada2f1 = xmg.dV./cada1f7;
cada1td1dV = cada2f1dV;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dV = -xmg.dV;
cada2f2 = uminus(xmg.f);
cada2f3dV = 2.*cada1f7.^(2-1).*cada1f7dV;
cada2f3 = cada1f7.^2;
cada2td1 = zeros(size(cada2f2dV,1),2);
cada2td1(:,2) = cada2f2dV./cada2f3;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = cada1f7dV(:,Gator2Data.Index27);
cada2td1 = cada2tf1.*cada2f4dV;
cada2td1(:,1) = cada2td1(:,1) + cada2f4.*cada1f7dVdV;
cada2f5dV = cada2td1;
cada2f5 = cada2f4.*cada1f7dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index28) = cada2f6dV;
cada2td1(:,2) = cada1td1dV(:,1);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f6;
term2.dVdV = cada1td1dV; term2.dV = cada1td1;
term2.f = xmg.f./cada1f7;
%User Line: term2 = xmg./(ones(size(t))+4./voverc+2./(voverc.^2));
cada2f1 = size(T.dV,1);
cada1td1 = zeros(cada2f1,3);
cada1td1(:,3) = T.dV;
cada2f1 = cada1td1(:,Gator1Data.Index9);
cada2f2dV = -D.dVdV;
cada2f2 = uminus(D.dV);
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada1td1dV = cada2f3dV;
cada1td1(:,Gator1Data.Index9) = cada2f3;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = T.f - D.f;
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,4);
cada1td1dV = cada1f1dVdV;
cada1td1(:,Gator1Data.Index10) = cada1f1dV;
cada2f1 = cada1td1(:,3);
cada2f2 = uminus(xmg.dV);
cada2f3 = cada2f1 + cada2f2;
cada1td1(:,3) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = cada1f1 - xmg.f;
cada2f1 = size(term1.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2tf1 = term1.dV(:,Gator2Data.Index29);
cada2td1 = cada2tf1.*term2.dV;
cada2td1(:,1) = cada2td1(:,1) + term2.f.*term1.dVdV;
cada2f1dV = cada2td1;
cada2f1 = term2.f.*term1.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada1tf1dV = term1.dV(:,Gator2Data.Index30);
cada1tf1 = term1.f(:,Gator1Data.Index11);
cada2td1 = zeros(size(cada1tf1dV,1),3);
cada2td1(:,Gator2Data.Index31) = term2.dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index32);
cada2td1 = cada2td1 + cada2tf1.*term2.dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*term2.dV;
cada2td1 = zeros(size(cada1td1dV,1),3);
cada2td1(:,Gator2Data.Index33) = cada1td1dV;
cada2td1 = cada2td1 + cada2f1dV;
cada1td1dV = cada2td1;
cada1td1 = cada1td1 + cada2f1;
cada1f3dVdV = cada1td1dV; cada1f3dV = cada1td1;
cada1f3 = term1.f.*term2.f;
cada1td1dV = cada1f2dVdV; cada1td1 = cada1f2dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index34);
cada2f1 = cada1td1(:,Gator1Data.Index12);
cada2f2dV = -cada1f3dVdV;
cada2f2 = uminus(cada1f3dV);
cada2td1 = zeros(size(cada2f1dV,1),4);
cada2td1(:,Gator2Data.Index35) = cada2f1dV;
cada2td1(:,Gator2Data.Index36) = cada2td1(:,Gator2Data.Index36) + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),6);
cada2td1(:,Gator2Data.Index37) = cada2f3dV;
cada2td1(:,Gator2Data.Index38) = cada1td1dV(:,Gator2Data.Index39);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index12) = cada2f3;
cada1f4dVdV = cada1td1dV; cada1f4dV = cada1td1;
cada1f4 = cada1f2 - cada1f3;
output(2).path.dVdV = cada1f4dVdV;
output(2).path.dV = cada1f4dV;
output(2).path.f = cada1f4;
%User Line: output(2).path = T-D-xmg-term1.*term2;
t.dV = input.phase(3).time.dV;
t.f = input.phase(3).time.f;
%User Line: t = input.phase(3).time;
x.dV = input.phase(3).state.dV;
x.f = input.phase(3).state.f;
%User Line: x = input.phase(3).state;
T.dV = input.phase(3).control.dV;
T.f = input.phase(3).control.f;
%User Line: T = input.phase(3).control;
h.dV = x.dV(:,1);
h.f = x.f(:,1);
%User Line: h = x(:,1);
v.dV = x.dV(:,2);
v.f = x.f(:,2);
%User Line: v = x(:,2);
m.dV = x.dV(:,3);
m.f = x.f(:,3);
%User Line: m = x(:,3);
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f1dVdV = v.dV.*cada2f2dV;
cada1f1dV = cada2f2.*v.dV;
cada1f1 = v.f.^2;
cada1f2dVdV = auxdata.dragk.*cada1f1dVdV;
cada1f2dV = auxdata.dragk*cada1f1dV;
cada1f2 = auxdata.dragk*cada1f1;
cada1f3dV = uminus(h.dV);
cada1f3 = uminus(h.f);
cada1f4dV = cada1f3dV/auxdata.H;
cada1f4 = cada1f3/auxdata.H;
cada2f1dV = exp(cada1f4).*cada1f4dV;
cada2f1 = exp(cada1f4);
cada1f5dVdV = cada1f4dV.*cada2f1dV;
cada1f5dV = cada2f1.*cada1f4dV;
cada1f5 = exp(cada1f4);
cada2f1 = size(cada1f2dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f5dV,1),2);
cada2td1(:,1) = cada1f2dV.*cada1f5dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f5.*cada1f2dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f5.*cada1f2dV;
cada1td1dV = cada2f1dV;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2td1 = zeros(size(cada1f2dV,1),2);
cada2td1(:,2) = cada1f5dV.*cada1f2dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f2.*cada1f5dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f2.*cada1f5dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index40) = cada2f3dV;
cada2td1(:,Gator2Data.Index41) = cada1td1dV(:,Gator2Data.Index42);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f3;
D.dVdV = cada1td1dV; D.dV = cada1td1;
D.f = cada1f2.*cada1f5;
%User Line: D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
hdot.dV = v.dV;
hdot.f = v.f;
%User Line: hdot = v;
cada2f1 = size(T.dV,1);
cada1td1 = zeros(cada2f1,3);
cada1td1(:,3) = T.dV;
cada2f1 = cada1td1(:,Gator1Data.Index13);
cada2f2dV = -D.dVdV;
cada2f2 = uminus(D.dV);
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada1td1dV = cada2f3dV;
cada1td1(:,Gator1Data.Index13) = cada2f3;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = T.f - D.f;
cada1tf1dV = m.dV(:,Gator2Data.Index43);
cada1tf1 = m.f(:,Gator1Data.Index14);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,4);
cada2tf1 = cada1tf1(:,Gator2Data.Index44);
cada2td1 = zeros(size(cada1f1dVdV,1),7);
cada2td1(:,Gator2Data.Index45) = cada1f1dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index46) = cada2td1(:,Gator2Data.Index46) + -cada1f1dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f1dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index15) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2f2dV = -cada1f1dV;
cada2f2 = uminus(cada1f1);
cada2f3dV = 2.*m.f.^(2-1).*m.dV;
cada2f3 = m.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index47);
cada2td1 = zeros(size(cada2f2dV,1),4);
cada2td1(:,Gator2Data.Index48) = cada2f2dV./cada2tf1;
cada2td1(:,3) = cada2td1(:,3) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = m.dV(:,Gator2Data.Index49);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*m.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),11);
cada2td1(:,Gator2Data.Index50) = cada2f6dV;
cada2td1(:,Gator2Data.Index51) = cada1td1dV(:,Gator2Data.Index52);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f6;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = cada1f1./m.f;
cada1f3 = size(t.f);
cada1f4 = ones(cada1f3);
cada1f5 = auxdata.g0*cada1f4;
vdot.dVdV = cada1f2dVdV; vdot.dV = cada1f2dV;
vdot.f = cada1f2 - cada1f5;
%User Line: vdot = (T-D)./m-auxdata.g0*ones(size(t));
cada1f1dV = uminus(T.dV);
cada1f1 = uminus(T.f);
mdot.dV = cada1f1dV/auxdata.c;
mdot.f = cada1f1/auxdata.c;
%User Line: mdot = -T./auxdata.c;
cada2f1 = size(hdot.f,1);
cada1td1 = zeros(cada2f1,6);
cada1td1(:,2) = hdot.dV;
cada1td1dV = vdot.dVdV;
cada1td1(:,Gator1Data.Index16) = vdot.dV;
cada1td1(:,6) = mdot.dV;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = [hdot.f vdot.f mdot.f];
output(3).dynamics.dVdV = cada1f1dVdV;
output(3).dynamics.dV = cada1f1dV;
output(3).dynamics.f = cada1f1;
%User Line: output(3).dynamics = [hdot, vdot, mdot];
%User Line: %------------------------------------%
%User Line: % End File:  goddardRocketContinuous %
%User Line: %------------------------------------%
cada2f1 = [3 5];
output(1).dynamics.dV_size = cada2f1;
output(1).dynamics.dV_location = Gator1Data.Index17;
cada2f1 = [3 5];
output(2).dynamics.dV_size = cada2f1;
output(2).dynamics.dV_location = Gator1Data.Index18;
cada2f1 = 5;
output(2).path.dV_size = cada2f1;
output(2).path.dV_location = Gator1Data.Index19;
cada2f1 = [3 5];
output(3).dynamics.dV_size = cada2f1;
output(3).dynamics.dV_location = Gator1Data.Index20;
output(1).dynamics.dVdV_size = [output(1).dynamics.dV_size,5];
output(1).dynamics.dVdV_location = [output(1).dynamics.dV_location(Gator2Data.Index53,:), Gator2Data.Index54];
output(2).dynamics.dVdV_size = [output(2).dynamics.dV_size,5];
output(2).dynamics.dVdV_location = [output(2).dynamics.dV_location(Gator2Data.Index55,:), Gator2Data.Index56];
output(2).path.dVdV_size = [output(2).path.dV_size,5];
output(2).path.dVdV_location = [output(2).path.dV_location(Gator2Data.Index57,:), Gator2Data.Index58];
output(3).dynamics.dVdV_size = [output(3).dynamics.dV_size,5];
output(3).dynamics.dVdV_location = [output(3).dynamics.dV_location(Gator2Data.Index59,:), Gator2Data.Index60];
end


function ADiGator_LoadData()
global ADiGator_goddardRocketContinuousADiGatorHes
ADiGator_goddardRocketContinuousADiGatorHes = load('goddardRocketContinuousADiGatorHes.mat');
return
end