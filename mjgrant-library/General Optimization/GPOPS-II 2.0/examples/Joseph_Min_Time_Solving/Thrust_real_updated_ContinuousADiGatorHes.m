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

function output = Thrust_real_updated_ContinuousADiGatorHes(input)
global ADiGator_Thrust_real_updated_ContinuousADiGatorHes
if isempty(ADiGator_Thrust_real_updated_ContinuousADiGatorHes); ADiGator_LoadData(); end
Gator1Data = ADiGator_Thrust_real_updated_ContinuousADiGatorHes.Thrust_real_updated_ContinuousADiGatorHes.Gator1Data;
Gator2Data = ADiGator_Thrust_real_updated_ContinuousADiGatorHes.Thrust_real_updated_ContinuousADiGatorHes.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: %--------------------------------------%
%User Line: % Begin File:  HMME_phases_Continuous %
%User Line: %--------------------------------------%
%User Line: % auxdata = input.auxdata;
cd0.f = input.auxdata.cd(1);
%User Line: cd0      = input.auxdata.cd(1);
cd1.f = input.auxdata.cd(2);
%User Line: cd1      = input.auxdata.cd(2);
cd2.f = input.auxdata.cd(3);
%User Line: cd2      = input.auxdata.cd(3);
cl0.f = input.auxdata.cl(1);
%User Line: cl0      = input.auxdata.cl(1);
cl1.f = input.auxdata.cl(2);
%User Line: cl1      = input.auxdata.cl(2);
%User Line: % keyboard
mu = input.auxdata.mu;
%User Line: mu  = input.auxdata.mu;
rho0 = input.auxdata.rho0;
%User Line: rho0 = input.auxdata.rho0;
H = input.auxdata.H;
%User Line: H = input.auxdata.H;
S = input.auxdata.S;
%User Line: S = input.auxdata.S;
%User Line: % Isp = input.auxdata.Isp;
g0 = input.auxdata.g0;
%User Line: g0 = input.auxdata.g0;
%User Line: % New parameters
y = input.auxdata.y;
%User Line: y  = input.auxdata.y;
R = input.auxdata.R;
%User Line: R  = input.auxdata.R;
T0 = input.auxdata.T0;
%User Line: T0  = input.auxdata.T0;
T4 = input.auxdata.T4;
%User Line: T4  = input.auxdata.T4;
cp = input.auxdata.cp;
%User Line: cp  = input.auxdata.cp;
hpr = input.auxdata.hpr;
%User Line: hpr  = input.auxdata.hpr;
Mc = input.auxdata.Mc;
%User Line: Mc = input.auxdata.Mc;
M1 = input.auxdata.M1;
%User Line: M1 = input.auxdata.M1;
M2 = input.auxdata.M2;
%User Line: M2 = input.auxdata.M2;
q1 = input.auxdata.q1;
%User Line: q1 = input.auxdata.q1;
q2 = input.auxdata.q2;
%User Line: q2 = input.auxdata.q2;
qscale = input.auxdata.qscale;
%User Line: qscale = input.auxdata.qscale;
T1 = input.auxdata.T1;
%User Line: T1 = input.auxdata.T1;
t.dV = input.phase.time.dV;
t.f = input.phase.time.f;
%User Line: t = input.phase.time;
x.dV = input.phase.state.dV;
x.f = input.phase.state.f;
%User Line: x = input.phase.state;
h.dV = x.dV(:,1);
h.f = x.f(:,1);
%User Line: h = x(:,1);
thetta.dV = x.dV(:,2);
thetta.f = x.f(:,2);
%User Line: thetta = x(:,2);
v.dV = x.dV(:,3);
v.f = x.f(:,3);
%User Line: v = x(:,3);
gam.dV = x.dV(:,4);
gam.f = x.f(:,4);
%User Line: gam = x(:,4);
mass.dV = x.dV(:,5);
mass.f = x.f(:,5);
%User Line: mass = x(:,5);
aoa.dV = x.dV(:,6);
aoa.f = x.f(:,6);
%User Line: aoa = x(:,6);
alfadot.dV = input.phase.control.dV(:,1);
alfadot.f = input.phase.control.f(:,1);
%User Line: alfadot  = input.phase.control(:,1);
A.dV = input.phase.control.dV(:,2);
A.f = input.phase.control.f(:,2);
%User Line: A = input.phase.control(:,2);
Re = input.auxdata.Re;
%User Line: Re = input.auxdata.Re;
%User Line: % Calculate stuff for phase 1
rad.dV = h.dV;
rad.f = h.f + Re;
%User Line: rad = h + Re;
cada1f1dV = uminus(h.dV);
cada1f1 = uminus(h.f);
cada1f2dV = cada1f1dV/H;
cada1f2 = cada1f1/H;
cada2f1dV = exp(cada1f2).*cada1f2dV;
cada2f1 = exp(cada1f2);
cada1f3dVdV = cada1f2dV.*cada2f1dV;
cada1f3dV = cada2f1.*cada1f2dV;
cada1f3 = exp(cada1f2);
rho.dVdV = rho0.*cada1f3dVdV;
rho.dV = rho0*cada1f3dV;
rho.f = rho0*cada1f3;
%User Line: rho = rho0*exp(-h/H);
cada1f1dV = cd1.f*aoa.dV;
cada1f1 = cd1.f*aoa.f;
cada1f2dV = cada1f1dV;
cada1f2 = cd0.f + cada1f1;
cada2f1dV = 1.*aoa.f.^(1-1).*aoa.dV;
cada2f1 = aoa.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f3dVdV = aoa.dV.*cada2f2dV;
cada1f3dV = cada2f2.*aoa.dV;
cada1f3 = aoa.f.^2;
cada1f4dVdV = cd2.f.*cada1f3dVdV;
cada1f4dV = cd2.f*cada1f3dV;
cada1f4 = cd2.f*cada1f3;
cada1td1 = cada1f2dV;
cada1td1dV = cada1f4dVdV;
cada1td1 = cada1td1 + cada1f4dV;
CD.dVdV = cada1td1dV; CD.dV = cada1td1;
CD.f = cada1f2 + cada1f4;
%User Line: CD       = cd0+cd1*aoa+cd2*aoa.^2;
cada1f1dV = cl1.f*aoa.dV;
cada1f1 = cl1.f*aoa.f;
CL.dV = cada1f1dV;
CL.f = cl0.f + cada1f1;
%User Line: CL       = cl0+cl1*aoa;
cada1f1dVdV = 0.5.*rho.dVdV;
cada1f1dV = 0.5*rho.dV;
cada1f1 = 0.5*rho.f;
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f2dVdV = v.dV.*cada2f2dV;
cada1f2dV = cada2f2.*v.dV;
cada1f2 = v.f.^2;
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f2dV,1),2);
cada2td1(:,2) = cada1f1dV.*cada1f2dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f2.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f2.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(cada1f1dV,1),2);
cada2td1(:,1) = cada1f2dV.*cada1f1dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f1.*cada1f2dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f1.*cada1f2dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index1) = cada2f3dV;
cada2td1(:,Gator2Data.Index2) = cada1td1dV(:,Gator2Data.Index3);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
dynamic_pressure.dVdV = cada1td1dV; dynamic_pressure.dV = cada1td1;
dynamic_pressure.f = cada1f1.*cada1f2;
%User Line: dynamic_pressure = 0.5*rho.*v.^2;
cada1f1dVdV = 0.5.*rho.dVdV;
cada1f1dV = 0.5*rho.dV;
cada1f1 = 0.5*rho.f;
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f2dVdV = v.dV.*cada2f2dV;
cada1f2dV = cada2f2.*v.dV;
cada1f2 = v.f.^2;
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f2dV,1),2);
cada2td1(:,2) = cada1f1dV.*cada1f2dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f2.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f2.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(cada1f1dV,1),2);
cada2td1(:,1) = cada1f2dV.*cada1f1dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f1.*cada1f2dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f1.*cada1f2dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index4) = cada2f3dV;
cada2td1(:,Gator2Data.Index5) = cada1td1dV(:,Gator2Data.Index6);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
q.dVdV = cada1td1dV; q.dV = cada1td1;
q.f = cada1f1.*cada1f2;
%User Line: q = 0.5*rho.*v.^2;
cada1f1dVdV = S.*dynamic_pressure.dVdV;
cada1f1dV = S*dynamic_pressure.dV;
cada1f1 = dynamic_pressure.f*S;
cada1tf1dV = CD.dV(:,Gator2Data.Index7);
cada1tf1 = CD.f(:,Gator1Data.Index1);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),6);
cada2td1(:,Gator2Data.Index8) = cada1f1dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index9);
cada2td1(:,Gator2Data.Index10) = cada2td1(:,Gator2Data.Index10) + cada2tf1.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index2) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = CD.dV(:,Gator2Data.Index11);
cada2td1 = zeros(size(cada1f1dV,1),3);
cada2td1(:,Gator2Data.Index12) = cada2tf1.*cada1f1dV;
cada2td1(:,3) = cada2td1(:,3) + cada1f1.*CD.dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f1.*CD.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),9);
cada2td1(:,Gator2Data.Index13) = cada2f3dV;
cada2td1(:,Gator2Data.Index14) = cada1td1dV(:,Gator2Data.Index15);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
D.dVdV = cada1td1dV; D.dV = cada1td1;
D.f = cada1f1.*CD.f;
%User Line: D = dynamic_pressure.*S.*CD;
cada1f1dVdV = S.*dynamic_pressure.dVdV;
cada1f1dV = S*dynamic_pressure.dV;
cada1f1 = dynamic_pressure.f*S;
cada1tf1dV = CL.dV(:,Gator2Data.Index16);
cada1tf1 = CL.f(:,Gator1Data.Index3);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),6);
cada2td1(:,Gator2Data.Index17) = cada1f1dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index18);
cada2td1(:,Gator2Data.Index19) = cada2td1(:,Gator2Data.Index19) + cada2tf1.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index4) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = CL.dV(:,Gator2Data.Index20);
cada2f2dV = cada2tf1.*cada1f1dV;
cada2f2 = cada1f1.*CL.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),8);
cada2td1(:,Gator2Data.Index21) = cada2f3dV;
cada2td1(:,Gator2Data.Index22) = cada1td1dV(:,Gator2Data.Index23);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
L.dVdV = cada1td1dV; L.dV = cada1td1;
L.f = cada1f1.*CL.f;
%User Line: L = dynamic_pressure.*S.*CL;
cada2f1dV = 1.*rad.f.^(1-1).*rad.dV;
cada2f1 = rad.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f1dVdV = rad.dV.*cada2f2dV;
cada1f1dV = cada2f2.*rad.dV;
cada1f1 = rad.f.^2;
cada2f1 = uminus(mu);
cada2f2dV = 2.*cada1f1.^(2-1).*cada1f1dV;
cada2f2 = cada1f1.^2;
cada2f3dV = -cada2f1./cada2f2.^2.*cada2f2dV;
cada2f3 = cada2f1./cada2f2;
cada2td1 = cada1f1dV.*cada2f3dV;
cada2td1 = cada2td1 + cada2f3.*cada1f1dVdV;
gravity.dVdV = cada2td1;
gravity.dV = cada2f3.*cada1f1dV;
gravity.f = mu./cada1f1;
%User Line: gravity = mu./rad.^2;
bank.f = 0;
%User Line: bank = 0;
cada2f1dV = -sin(gam.f).*gam.dV;
cada2f1 = cos(gam.f);
sgam.dVdV = gam.dV.*cada2f1dV;
sgam.dV = cada2f1.*gam.dV;
sgam.f = sin(gam.f);
%User Line: sgam = sin(gam);
cada2f1dV = cos(gam.f).*gam.dV;
cada2f1 = sin(gam.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cgam.dVdV = gam.dV.*cada2f2dV;
cgam.dV = cada2f2.*gam.dV;
cgam.f = cos(gam.f);
%User Line: cgam = cos(gam);
cbank.f = cos(bank.f);
%User Line: cbank = cos(bank);
%User Line: % New equations needed
cada1f1 = y*R;
cada1f2 = cada1f1*T0;
va.f = sqrt(cada1f2);
%User Line: va = sqrt(y*R*T0);
M.dV = v.dV/va.f;
M.f = v.f/va.f;
%User Line: M = v./va;
cada1f1dV = uminus(h.dV);
cada1f1 = uminus(h.f);
cada1f2dV = cada1f1dV/H;
cada1f2 = cada1f1/H;
cada2f1dV = exp(cada1f2).*cada1f2dV;
cada2f1 = exp(cada1f2);
cada1f3dVdV = cada1f2dV.*cada2f1dV;
cada1f3dV = cada2f1.*cada1f2dV;
cada1f3 = exp(cada1f2);
rho.dVdV = rho0.*cada1f3dVdV;
rho.dV = rho0*cada1f3dV;
rho.f = rho0*cada1f3;
%User Line: rho = rho0*exp(-h/H);
cada2f1 = size(rho.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(A.dV,1),2);
cada2td1(:,2) = rho.dV.*A.dV;
cada2td1(:,1) = cada2td1(:,1) + A.f.*rho.dVdV;
cada2f1dV = cada2td1;
cada2f1 = A.f.*rho.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2f2dV = A.dV.*rho.dV;
cada2f2 = rho.f.*A.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,2) = cada2f3dV;
cada2td1(:,Gator2Data.Index24) = cada1td1dV(:,Gator2Data.Index25);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = rho.f.*A.f;
cada1tf1dV = v.dV(:,Gator2Data.Index26);
cada1tf1 = v.f(:,Gator1Data.Index5);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),5);
cada2td1(:,Gator2Data.Index27) = cada1f1dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index28);
cada2td1(:,Gator2Data.Index29) = cada2td1(:,Gator2Data.Index29) + cada2tf1.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index6) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2tf1 = v.dV(:,Gator2Data.Index30);
cada2f2dV = cada2tf1.*cada1f1dV;
cada2f2 = cada1f1.*v.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),7);
cada2td1(:,Gator2Data.Index31) = cada2f3dV;
cada2td1(:,Gator2Data.Index32) = cada1td1dV(:,Gator2Data.Index33);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
madot.dVdV = cada1td1dV; madot.dV = cada1td1;
madot.f = cada1f1.*v.f;
%User Line: madot = rho.*A.*v;
cada1f1 = y - 1;
cada1f2 = Mc^2;
cada1f3 = cada1f1*cada1f2;
cada1f4 = cada1f3/2;
cada1f5 = 1 + cada1f4;
cada1f6 = T4*cada1f5;
Tau_lam.f = cada1f6/T0;
%User Line: Tau_lam = T4*(1+(y-1)*Mc^2/2)/T0;
cada1f1 = y - 1;
cada2f1dV = 1.*M.f.^(1-1).*M.dV;
cada2f1 = M.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f2dVdV = M.dV.*cada2f2dV;
cada1f2dV = cada2f2.*M.dV;
cada1f2 = M.f.^2;
cada1f3dVdV = cada1f1.*cada1f2dVdV;
cada1f3dV = cada1f1*cada1f2dV;
cada1f3 = cada1f1*cada1f2;
cada1f4dVdV = cada1f3dVdV./2;
cada1f4dV = cada1f3dV/2;
cada1f4 = cada1f3/2;
Tau_r.dVdV = cada1f4dVdV; Tau_r.dV = cada1f4dV;
Tau_r.f = 1 + cada1f4;
%User Line: Tau_r = (1+(y-1).*M.^2/2);
%User Line: %Tau_lam = Tau_r +(Tau_lam - Tau_r).*(atan(qscale*(q-q1))-atan(qscale*(q-q2)))/pi;
cada1f1 = cp*T0;
cada1f2 = cada1f1/hpr;
cada1f3dVdV = -Tau_r.dVdV;
cada1f3dV = uminus(Tau_r.dV);
cada1f3 = Tau_lam.f - Tau_r.f;
f.dVdV = cada1f2.*cada1f3dVdV;
f.dV = cada1f2*cada1f3dV;
f.f = cada1f2*cada1f3;
%User Line: f = cp*T0/hpr*(Tau_lam-Tau_r);
cada2f1 = uminus(Tau_lam.f);
cada2f2dV = 2.*Tau_r.f.^(2-1).*Tau_r.dV;
cada2f2 = Tau_r.f.^2;
cada2f3dV = -cada2f1./cada2f2.^2.*cada2f2dV;
cada2f3 = cada2f1./cada2f2;
cada2td1 = Tau_r.dV.*cada2f3dV;
cada2td1 = cada2td1 + cada2f3.*Tau_r.dVdV;
cada1f1dVdV = cada2td1;
cada1f1dV = cada2f3.*Tau_r.dV;
cada1f1 = Tau_lam.f./Tau_r.f;
cada2f1dV = (1/2)./sqrt(cada1f1).*cada1f1dV;
cada2f1dV(cada1f1 == 0 & cada1f1dV == 0) = 0;
cada2f1 = sqrt(cada1f1);
cada2f2dV = -0.5./cada2f1.^2.*cada2f1dV;
cada2f2 = 0.5./cada2f1;
cada2td1 = cada1f1dV.*cada2f2dV;
cada2td1 = cada2td1 + cada2f2.*cada1f1dVdV;
cada1f2dVdV = cada2td1;
cada1f2dV = cada2f2.*cada1f1dV;
cada2f1 = eq(cada1f1,0);
cada2f2 = eq(cada1f1dV,0);
cada2f3 = and(cada2f1,cada2f2);
cada2td2 = cada1f2dVdV;
cada2tind1 = cada2f3(:,1);
cada2td2(cada2tind1) = 0;
cada1f2dVdV = cada2td2;
cada1f2dV(cada2f3) = 0;
cada1f2 = sqrt(cada1f1);
cada1f3dVdV = cada1f2dVdV; cada1f3dV = cada1f2dV;
cada1f3 = cada1f2 - 1;
cada1td1dV = v.dV.*cada1f3dV;
cada1td1 = cada1f3.*v.dV;
cada2td1 = cada1f3dV.*v.dV;
cada2td1 = cada2td1 + v.f.*cada1f3dVdV;
cada2f1dV = cada2td1;
cada2f1 = v.f.*cada1f3dV;
cada2td1 = cada1td1dV;
cada2td1 = cada2td1 + cada2f1dV;
cada1td1dV = cada2td1;
cada1td1 = cada1td1 + cada2f1;
ST.dVdV = cada1td1dV; ST.dV = cada1td1;
ST.f = v.f.*cada1f3;
%User Line: ST = v.*(sqrt(Tau_lam./Tau_r)-1);
cada1f1dVdV = g0.*f.dVdV;
cada1f1dV = g0*f.dV;
cada1f1 = f.f*g0;
cada2td1 = ST.dVdV./cada1f1;
cada2td1 = cada2td1 + -ST.dV./cada1f1.^2.*cada1f1dV;
cada1td1dV = cada2td1;
cada1td1 = ST.dV./cada1f1;
cada2f1dV = -ST.dV;
cada2f1 = uminus(ST.f);
cada2f2dV = 2.*cada1f1.^(2-1).*cada1f1dV;
cada2f2 = cada1f1.^2;
cada2td1 = cada2f1dV./cada2f2;
cada2td1 = cada2td1 + -cada2f1./cada2f2.^2.*cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1./cada2f2;
cada2td1 = cada1f1dV.*cada2f3dV;
cada2td1 = cada2td1 + cada2f3.*cada1f1dVdV;
cada2f4dV = cada2td1;
cada2f4 = cada2f3.*cada1f1dV;
cada2td1 = cada1td1dV;
cada2td1 = cada2td1 + cada2f4dV;
cada1td1dV = cada2td1;
cada1td1 = cada1td1 + cada2f4;
Isp.dVdV = cada1td1dV; Isp.dV = cada1td1;
Isp.f = ST.f./cada1f1;
%User Line: Isp = ST./(f*g0);
%User Line: % % keyboard
%User Line: % for i = 1:1:length(ST)
%User Line: %     % M_info = M(i);
%User Line: %     q_info = dynamic_pressure(i);
%User Line: %     if (q_info<40000 || q_info>150000) % M_info<4 || M_info>9.75 ||
%User Line: %         T(i) = 0;
%User Line: %     else
%User Line: %       T(i) = madot(i)*ST(i);
%User Line: %     end
%User Line: % end
%User Line: % T = T';
cada1tf1dV = ST.dV(:,Gator2Data.Index34);
cada1tf1 = ST.f(:,Gator1Data.Index7);
cada2td1 = zeros(size(cada1tf1dV,1),8);
cada2td1(:,Gator2Data.Index35) = madot.dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index36);
cada2td1(:,Gator2Data.Index37) = cada2td1(:,Gator2Data.Index37) + cada2tf1.*madot.dVdV;
cada1td1dV = cada2td1;
cada1td1 = cada1tf1.*madot.dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index38);
cada2f1 = cada1td1(:,2);
cada2tf1 = ST.dV(:,Gator2Data.Index39);
cada2td1 = cada2tf1.*madot.dV;
cada2td1(:,2) = cada2td1(:,2) + madot.f.*ST.dVdV;
cada2f2dV = cada2td1;
cada2f2 = madot.f.*ST.dV;
cada2td1 = cada2f1dV;
cada2td1 = cada2td1 + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),8);
cada2td1(:,Gator2Data.Index40) = cada2f3dV;
cada2td1(:,Gator2Data.Index41) = cada1td1dV(:,Gator2Data.Index42);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
T.dVdV = cada1td1dV; T.dV = cada1td1;
T.f = madot.f.*ST.f;
%User Line: T = madot.*ST;
%User Line: % T = T.*(atan(1e17*(M-M1))+atan(1e17*(M-M2)))/pi;
%User Line: % Tmax = input.auxdata.Tmax;
%User Line: % T = Tmax*A;
%User Line: % path = dynamic_pressure;
%User Line: % output.path = path;
cada2f1 = size(v.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dV = v.dV.*sgam.dV;
cada2f1 = sgam.f.*v.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dV,1),2);
cada2td1(:,1) = sgam.dV.*v.dV;
cada2td1(:,2) = cada2td1(:,2) + v.f.*sgam.dVdV;
cada2f2dV = cada2td1;
cada2f2 = v.f.*sgam.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index43) = cada2f3dV;
cada2td1(:,2) = cada1td1dV(:,1);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
hdot.dVdV = cada1td1dV; hdot.dV = cada1td1;
hdot.f = v.f.*sgam.f;
%User Line: hdot   = v.*sgam;
cada2f1 = size(v.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dV = v.dV.*cgam.dV;
cada2f1 = cgam.f.*v.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dV,1),2);
cada2td1(:,1) = cgam.dV.*v.dV;
cada2td1(:,2) = cada2td1(:,2) + v.f.*cgam.dVdV;
cada2f2dV = cada2td1;
cada2f2 = v.f.*cgam.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index44) = cada2f3dV;
cada2td1(:,2) = cada1td1dV(:,1);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f1dVdV = cada1td1dV; cada1f1dV = cada1td1;
cada1f1 = v.f.*cgam.f;
cada1tf1dV = rad.dV(:,Gator2Data.Index45);
cada1tf1 = rad.f(:,Gator1Data.Index8);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,3);
cada2tf1 = cada1tf1(:,Gator2Data.Index46);
cada2td1 = zeros(size(cada1f1dVdV,1),5);
cada2td1(:,Gator2Data.Index47) = cada1f1dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index48) = cada2td1(:,Gator2Data.Index48) + -cada1f1dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f1dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index9) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dV = -cada1f1dV;
cada2f2 = uminus(cada1f1);
cada2f3dV = 2.*rad.f.^(2-1).*rad.dV;
cada2f3 = rad.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index49);
cada2td1 = zeros(size(cada2f2dV,1),3);
cada2td1(:,Gator2Data.Index50) = cada2f2dV./cada2tf1;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = rad.dV(:,Gator2Data.Index51);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*rad.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),8);
cada2td1(:,Gator2Data.Index52) = cada2f6dV;
cada2td1(:,Gator2Data.Index53) = cada1td1dV(:,Gator2Data.Index54);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f6;
thettadot.dVdV = cada1td1dV; thettadot.dV = cada1td1;
thettadot.f = cada1f1./rad.f;
%User Line: thettadot = v.*cgam./rad;
cada2f1dV = cos(aoa.f).*aoa.dV;
cada2f1 = sin(aoa.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f1dVdV = aoa.dV.*cada2f2dV;
cada1f1dV = cada2f2.*aoa.dV;
cada1f1 = cos(aoa.f);
cada1tf1dV = cada1f1dV(:,Gator2Data.Index55);
cada1tf1 = cada1f1(:,Gator1Data.Index10);
cada2f1 = size(T.dV,1);
cada1td1 = zeros(cada2f1,4);
cada2td1 = zeros(size(cada1tf1dV,1),11);
cada2td1(:,Gator2Data.Index56) = T.dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index57);
cada2td1(:,Gator2Data.Index58) = cada2td1(:,Gator2Data.Index58) + cada2tf1.*T.dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*T.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index11) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = cada1f1dV(:,Gator2Data.Index59);
cada2td1 = zeros(size(T.dV,1),4);
cada2td1(:,Gator2Data.Index60) = cada2tf1.*T.dV;
cada2td1(:,3) = cada2td1(:,3) + T.f.*cada1f1dVdV;
cada2f2dV = cada2td1;
cada2f2 = T.f.*cada1f1dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),15);
cada2td1(:,Gator2Data.Index61) = cada2f3dV;
cada2td1(:,Gator2Data.Index62) = cada1td1dV(:,Gator2Data.Index63);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = T.f.*cada1f1;
cada1td1dV = cada1f2dVdV; cada1td1 = cada1f2dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index64);
cada2f1 = cada1td1(:,Gator1Data.Index12);
cada2f2dV = -D.dVdV;
cada2f2 = uminus(D.dV);
cada2td1 = cada2f1dV;
cada2td1(:,Gator2Data.Index65) = cada2td1(:,Gator2Data.Index65) + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),15);
cada2td1(:,Gator2Data.Index66) = cada2f3dV;
cada2td1(:,Gator2Data.Index67) = cada1td1dV(:,Gator2Data.Index68);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index12) = cada2f3;
cada1f3dVdV = cada1td1dV; cada1f3dV = cada1td1;
cada1f3 = cada1f2 - D.f;
cada1tf1dV = mass.dV(:,Gator2Data.Index69);
cada1tf1 = mass.f(:,Gator1Data.Index13);
cada2f1 = size(cada1f3dV,1);
cada1td1 = zeros(cada2f1,5);
cada2tf1 = cada1tf1(:,Gator2Data.Index70);
cada2td1 = zeros(size(cada1f3dVdV,1),19);
cada2td1(:,Gator2Data.Index71) = cada1f3dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index72) = cada2td1(:,Gator2Data.Index72) + -cada1f3dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f3dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index14) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2f2dV = -cada1f3dV;
cada2f2 = uminus(cada1f3);
cada2f3dV = 2.*mass.f.^(2-1).*mass.dV;
cada2f3 = mass.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index73);
cada2td1 = zeros(size(cada2f2dV,1),5);
cada2td1(:,Gator2Data.Index74) = cada2f2dV./cada2tf1;
cada2td1(:,3) = cada2td1(:,3) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = mass.dV(:,Gator2Data.Index75);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*mass.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),24);
cada2td1(:,Gator2Data.Index76) = cada2f6dV;
cada2td1(:,Gator2Data.Index77) = cada1td1dV(:,Gator2Data.Index78);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f6;
cada1f4dVdV = cada1td1dV; cada1f4dV = cada1td1;
cada1f4 = cada1f3./mass.f;
cada2f1 = size(gravity.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(sgam.dV,1),2);
cada2td1(:,2) = gravity.dV.*sgam.dV;
cada2td1(:,1) = cada2td1(:,1) + sgam.f.*gravity.dVdV;
cada2f1dV = cada2td1;
cada2f1 = sgam.f.*gravity.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(gravity.dV,1),2);
cada2td1(:,1) = sgam.dV.*gravity.dV;
cada2td1(:,2) = cada2td1(:,2) + gravity.f.*sgam.dVdV;
cada2f2dV = cada2td1;
cada2f2 = gravity.f.*sgam.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index79) = cada2f3dV;
cada2td1(:,Gator2Data.Index80) = cada1td1dV(:,Gator2Data.Index81);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f5dVdV = cada1td1dV; cada1f5dV = cada1td1;
cada1f5 = gravity.f.*sgam.f;
cada2f1 = size(cada1f4dV,1);
cada1td1 = zeros(cada2f1,6);
cada1td1dV = cada1f4dVdV;
cada1td1(:,Gator1Data.Index15) = cada1f4dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index82);
cada2f1 = cada1td1(:,Gator1Data.Index16);
cada2f2dV = -cada1f5dVdV;
cada2f2 = uminus(cada1f5dV);
cada2td1 = zeros(size(cada2f1dV,1),8);
cada2td1(:,Gator2Data.Index83) = cada2f1dV;
cada2td1(:,Gator2Data.Index84) = cada2td1(:,Gator2Data.Index84) + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),27);
cada2td1(:,Gator2Data.Index85) = cada2f3dV;
cada2td1(:,Gator2Data.Index86) = cada1td1dV(:,Gator2Data.Index87);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index16) = cada2f3;
vdot.dVdV = cada1td1dV; vdot.dV = cada1td1;
vdot.f = cada1f4 - cada1f5;
%User Line: vdot = (T.*cos(aoa)-D)./mass-gravity.*sgam;
cada2f1dV = -sin(aoa.f).*aoa.dV;
cada2f1 = cos(aoa.f);
cada1f1dVdV = aoa.dV.*cada2f1dV;
cada1f1dV = cada2f1.*aoa.dV;
cada1f1 = sin(aoa.f);
cada1tf1dV = cada1f1dV(:,Gator2Data.Index88);
cada1tf1 = cada1f1(:,Gator1Data.Index17);
cada2f1 = size(T.dV,1);
cada1td1 = zeros(cada2f1,4);
cada2td1 = zeros(size(cada1tf1dV,1),11);
cada2td1(:,Gator2Data.Index89) = T.dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index90);
cada2td1(:,Gator2Data.Index91) = cada2td1(:,Gator2Data.Index91) + cada2tf1.*T.dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*T.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index18) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = cada1f1dV(:,Gator2Data.Index92);
cada2td1 = zeros(size(T.dV,1),4);
cada2td1(:,Gator2Data.Index93) = cada2tf1.*T.dV;
cada2td1(:,3) = cada2td1(:,3) + T.f.*cada1f1dVdV;
cada2f2dV = cada2td1;
cada2f2 = T.f.*cada1f1dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),15);
cada2td1(:,Gator2Data.Index94) = cada2f3dV;
cada2td1(:,Gator2Data.Index95) = cada1td1dV(:,Gator2Data.Index96);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = T.f.*cada1f1;
cada1f3dVdV = cbank.f.*L.dVdV;
cada1f3dV = cbank.f*L.dV;
cada1f3 = L.f*cbank.f;
cada1td1dV = cada1f2dVdV; cada1td1 = cada1f2dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index97);
cada2f1 = cada1td1(:,Gator1Data.Index19);
cada2td1 = cada2f1dV;
cada2td1(:,Gator2Data.Index98) = cada2td1(:,Gator2Data.Index98) + cada1f3dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada2f1 + cada1f3dV;
cada2td1 = zeros(size(cada1td1,1),15);
cada2td1(:,Gator2Data.Index99) = cada2f2dV;
cada2td1(:,Gator2Data.Index100) = cada1td1dV(:,Gator2Data.Index101);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index19) = cada2f2;
cada1f4dVdV = cada1td1dV; cada1f4dV = cada1td1;
cada1f4 = cada1f2 + cada1f3;
cada1tf1dV = mass.dV(:,Gator2Data.Index102);
cada1tf1 = mass.f(:,Gator1Data.Index20);
cada2f1 = size(cada1f4dV,1);
cada1td1 = zeros(cada2f1,5);
cada2tf1 = cada1tf1(:,Gator2Data.Index103);
cada2td1 = zeros(size(cada1f4dVdV,1),19);
cada2td1(:,Gator2Data.Index104) = cada1f4dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index105) = cada2td1(:,Gator2Data.Index105) + -cada1f4dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f4dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index21) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2f2dV = -cada1f4dV;
cada2f2 = uminus(cada1f4);
cada2f3dV = 2.*mass.f.^(2-1).*mass.dV;
cada2f3 = mass.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index106);
cada2td1 = zeros(size(cada2f2dV,1),5);
cada2td1(:,Gator2Data.Index107) = cada2f2dV./cada2tf1;
cada2td1(:,3) = cada2td1(:,3) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = mass.dV(:,Gator2Data.Index108);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*mass.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),24);
cada2td1(:,Gator2Data.Index109) = cada2f6dV;
cada2td1(:,Gator2Data.Index110) = cada1td1dV(:,Gator2Data.Index111);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f6;
cada1f5dVdV = cada1td1dV; cada1f5dV = cada1td1;
cada1f5 = cada1f4./mass.f;
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f6dVdV = v.dV.*cada2f2dV;
cada1f6dV = cada2f2.*v.dV;
cada1f6 = v.f.^2;
cada2f1 = size(cada1f6dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f6dVdV,1),2);
cada2td1(:,2) = cada1f6dVdV./rad.f;
cada2td1(:,1) = cada2td1(:,1) + -cada1f6dV./rad.f.^2.*rad.dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f6dV./rad.f;
cada1td1dV = cada2f1dV;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dV = -cada1f6dV;
cada2f2 = uminus(cada1f6);
cada2f3dV = 2.*rad.f.^(2-1).*rad.dV;
cada2f3 = rad.f.^2;
cada2td1 = zeros(size(cada2f2dV,1),2);
cada2td1(:,2) = cada2f2dV./cada2f3;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = rad.dV(:,Gator2Data.Index112);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*rad.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index113) = cada2f6dV;
cada2td1(:,Gator2Data.Index114) = cada1td1dV(:,Gator2Data.Index115);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f6;
cada1f7dVdV = cada1td1dV; cada1f7dV = cada1td1;
cada1f7 = cada1f6./rad.f;
cada2f1 = size(gravity.dV,1);
cada1td1 = zeros(cada2f1,2);
cada1td1dV = gravity.dVdV;
cada1td1(:,1) = gravity.dV;
cada2f1dV = -cada1f7dVdV;
cada2f1 = uminus(cada1f7dV);
cada2td1 = zeros(size(cada1td1dV,1),4);
cada2td1(:,1) = cada1td1dV;
cada2td1 = cada2td1 + cada2f1dV;
cada1td1dV = cada2td1;
cada1td1 = cada1td1 + cada2f1;
cada1f8dVdV = cada1td1dV; cada1f8dV = cada1td1;
cada1f8 = gravity.f - cada1f7;
cada2f1 = size(cgam.dV,1);
cada1td1 = zeros(cada2f1,3);
cada2tf1 = cgam.dV(:,Gator2Data.Index116);
cada2td1 = zeros(size(cada1f8dV,1),3);
cada2td1(:,Gator2Data.Index117) = cada2tf1.*cada1f8dV;
cada2td1(:,3) = cada2td1(:,3) + cada1f8.*cgam.dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f8.*cgam.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,3) = cada2f1;
cada1tf1dV = cgam.dV(:,Gator2Data.Index118);
cada1tf1 = cgam.f(:,Gator1Data.Index22);
cada2f1 = cada1td1(:,Gator1Data.Index23);
cada2td1 = zeros(size(cada1tf1dV,1),6);
cada2td1(:,Gator2Data.Index119) = cada1f8dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index120);
cada2td1(:,Gator2Data.Index121) = cada2td1(:,Gator2Data.Index121) + cada2tf1.*cada1f8dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1tf1.*cada1f8dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),9);
cada2td1(:,Gator2Data.Index122) = cada2f3dV;
cada2td1(:,Gator2Data.Index123) = cada1td1dV(:,Gator2Data.Index124);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index23) = cada2f3;
cada1f9dVdV = cada1td1dV; cada1f9dV = cada1td1;
cada1f9 = cgam.f.*cada1f8;
cada2f1 = size(cada1f5dV,1);
cada1td1 = zeros(cada2f1,6);
cada1td1dV = cada1f5dVdV;
cada1td1(:,Gator1Data.Index24) = cada1f5dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index125);
cada2f1 = cada1td1(:,Gator1Data.Index25);
cada2f2dV = -cada1f9dVdV;
cada2f2 = uminus(cada1f9dV);
cada2td1 = zeros(size(cada2f1dV,1),15);
cada2td1(:,Gator2Data.Index126) = cada2f1dV;
cada2td1(:,Gator2Data.Index127) = cada2td1(:,Gator2Data.Index127) + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),29);
cada2td1(:,Gator2Data.Index128) = cada2f3dV;
cada2td1(:,Gator2Data.Index129) = cada1td1dV(:,Gator2Data.Index130);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index25) = cada2f3;
cada1f10dVdV = cada1td1dV; cada1f10dV = cada1td1;
cada1f10 = cada1f5 - cada1f9;
cada1tf1dV = v.dV(:,Gator2Data.Index131);
cada1tf1 = v.f(:,Gator1Data.Index26);
cada2tf1 = cada1tf1(:,Gator2Data.Index132);
cada2td1 = cada1f10dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index133) = cada2td1(:,Gator2Data.Index133) + -cada1f10dV./cada1tf1.^2.*cada1tf1dV;
cada1td1dV = cada2td1;
cada1td1 = cada1f10dV./cada1tf1;
cada2f1dV = cada1td1dV(:,Gator2Data.Index134);
cada2f1 = cada1td1(:,2);
cada2f2dV = -cada1f10dV;
cada2f2 = uminus(cada1f10);
cada2f3dV = 2.*v.f.^(2-1).*v.dV;
cada2f3 = v.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index135);
cada2td1 = cada2f2dV./cada2tf1;
cada2td1(:,2) = cada2td1(:,2) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = v.dV(:,Gator2Data.Index136);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*v.dV;
cada2td1 = cada2f1dV;
cada2td1 = cada2td1 + cada2f5dV;
cada2f6dV = cada2td1;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),29);
cada2td1(:,Gator2Data.Index137) = cada2f6dV;
cada2td1(:,Gator2Data.Index138) = cada1td1dV(:,Gator2Data.Index139);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f6;
gamdot.dVdV = cada1td1dV; gamdot.dV = cada1td1;
gamdot.f = cada1f10./v.f;
%User Line: gamdot   = ((T.*sin(aoa)+L.*cbank)./mass-cgam.*(gravity-v.^2./rad))./v;
cada1f1dVdV = -madot.dVdV;
cada1f1dV = uminus(madot.dV);
cada1f1 = uminus(madot.f);
cada1tf1dV = f.dV(:,Gator2Data.Index140);
cada1tf1 = f.f(:,Gator1Data.Index27);
cada2td1 = zeros(size(cada1tf1dV,1),8);
cada2td1(:,Gator2Data.Index141) = cada1f1dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index142);
cada2td1(:,Gator2Data.Index143) = cada2td1(:,Gator2Data.Index143) + cada2tf1.*cada1f1dVdV;
cada1td1dV = cada2td1;
cada1td1 = cada1tf1.*cada1f1dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index144);
cada2f1 = cada1td1(:,2);
cada2tf1 = f.dV(:,Gator2Data.Index145);
cada2td1 = cada2tf1.*cada1f1dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f1.*f.dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f1.*f.dV;
cada2td1 = cada2f1dV;
cada2td1 = cada2td1 + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),8);
cada2td1(:,Gator2Data.Index146) = cada2f3dV;
cada2td1(:,Gator2Data.Index147) = cada1td1dV(:,Gator2Data.Index148);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
massdot.dVdV = cada1td1dV; massdot.dV = cada1td1;
massdot.f = cada1f1.*f.f;
%User Line: massdot = -madot.*f;
aoadot.dV = alfadot.dV;
aoadot.f = alfadot.f;
%User Line: aoadot = alfadot;
cada2f1 = size(hdot.f,1);
cada1td1 = zeros(cada2f1,21);
cada1td1dV = hdot.dVdV;
cada1td1(:,Gator1Data.Index28) = hdot.dV;
cada2td1 = zeros(size(cada1td1,1),11);
cada2td1(:,Gator2Data.Index149) = thettadot.dVdV;
cada2td1(:,Gator2Data.Index150) = cada1td1dV(:,Gator2Data.Index151);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index29) = thettadot.dV;
cada2td1 = zeros(size(cada1td1,1),38);
cada2td1(:,Gator2Data.Index152) = vdot.dVdV;
cada2td1(:,Gator2Data.Index153) = cada1td1dV(:,Gator2Data.Index154);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index30) = vdot.dV;
cada2td1 = zeros(size(cada1td1,1),67);
cada2td1(:,Gator2Data.Index155) = gamdot.dVdV;
cada2td1(:,Gator2Data.Index156) = cada1td1dV(:,Gator2Data.Index157);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index31) = gamdot.dV;
cada2td1 = zeros(size(cada1td1,1),75);
cada2td1(:,Gator2Data.Index158) = massdot.dVdV;
cada2td1(:,Gator2Data.Index159) = cada1td1dV(:,Gator2Data.Index160);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index32) = massdot.dV;
cada1td1(:,18) = aoadot.dV;
output.dynamics.dVdV = cada1td1dV; output.dynamics.dV = cada1td1;
output.dynamics.f = [hdot.f thettadot.f vdot.f gamdot.f massdot.f aoadot.f];
%User Line: output.dynamics  = [hdot, thettadot, vdot, gamdot, massdot, aoadot];
%User Line: %------------------------------------%
%User Line: % End File:  goddardRocketContinuous %
%User Line: %------------------------------------%
output.dynamics.dV_size = [6 9];
output.dynamics.dV_location = Gator1Data.Index33;
output.dynamics.dVdV_size = [output.dynamics.dV_size,9];
output.dynamics.dVdV_location = [output.dynamics.dV_location(Gator2Data.Index161,:), Gator2Data.Index162];
end


function ADiGator_LoadData()
global ADiGator_Thrust_real_updated_ContinuousADiGatorHes
ADiGator_Thrust_real_updated_ContinuousADiGatorHes = load('Thrust_real_updated_ContinuousADiGatorHes.mat');
return
end