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

function output = Thrust_real_updated_ContinuousADiGatorGrd(input)
global ADiGator_Thrust_real_updated_ContinuousADiGatorGrd
if isempty(ADiGator_Thrust_real_updated_ContinuousADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_Thrust_real_updated_ContinuousADiGatorGrd.Thrust_real_updated_ContinuousADiGatorGrd.Gator1Data;
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
t.dV = input.phase.time.dV; t.f = input.phase.time.f;
%User Line: t = input.phase.time;
x.dV = input.phase.state.dV; x.f = input.phase.state.f;
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
cada1f1dV = -h.dV;
cada1f1 = uminus(h.f);
cada1f2dV = cada1f1dV./H;
cada1f2 = cada1f1/H;
cada1f3dV = exp(cada1f2).*cada1f2dV;
cada1f3 = exp(cada1f2);
rho.dV = rho0.*cada1f3dV;
rho.f = rho0*cada1f3;
%User Line: rho = rho0*exp(-h/H);
cada1f1dV = cd1.f.*aoa.dV;
cada1f1 = cd1.f*aoa.f;
cada1f2dV = cada1f1dV;
cada1f2 = cd0.f + cada1f1;
cada1f3dV = 2.*aoa.f.^(2-1).*aoa.dV;
cada1f3 = aoa.f.^2;
cada1f4dV = cd2.f.*cada1f3dV;
cada1f4 = cd2.f*cada1f3;
cada1td1 = cada1f2dV;
cada1td1 = cada1td1 + cada1f4dV;
CD.dV = cada1td1;
CD.f = cada1f2 + cada1f4;
%User Line: CD       = cd0+cd1*aoa+cd2*aoa.^2;
cada1f1dV = cl1.f.*aoa.dV;
cada1f1 = cl1.f*aoa.f;
CL.dV = cada1f1dV;
CL.f = cl0.f + cada1f1;
%User Line: CL       = cl0+cl1*aoa;
cada1f1dV = 0.5.*rho.dV;
cada1f1 = 0.5*rho.f;
cada1f2dV = 2.*v.f.^(2-1).*v.dV;
cada1f2 = v.f.^2;
cada1td1 = zeros(size(cada1f1dV,1),2);
cada1td1(:,1) = cada1f2.*cada1f1dV;
cada1td1(:,2) = cada1td1(:,2) + cada1f1.*cada1f2dV;
dynamic_pressure.dV = cada1td1;
dynamic_pressure.f = cada1f1.*cada1f2;
%User Line: dynamic_pressure = 0.5*rho.*v.^2;
cada1f1dV = 0.5.*rho.dV;
cada1f1 = 0.5*rho.f;
cada1f2dV = 2.*v.f.^(2-1).*v.dV;
cada1f2 = v.f.^2;
cada1td1 = zeros(size(cada1f1dV,1),2);
cada1td1(:,1) = cada1f2.*cada1f1dV;
cada1td1(:,2) = cada1td1(:,2) + cada1f1.*cada1f2dV;
q.dV = cada1td1;
q.f = cada1f1.*cada1f2;
%User Line: q = 0.5*rho.*v.^2;
cada1f1dV = S.*dynamic_pressure.dV;
cada1f1 = dynamic_pressure.f*S;
cada1tf1 = CD.f(:,Gator1Data.Index1);
cada1td1 = zeros(size(cada1f1dV,1),3);
cada1td1(:,Gator1Data.Index2) = cada1tf1.*cada1f1dV;
cada1td1(:,3) = cada1td1(:,3) + cada1f1.*CD.dV;
D.dV = cada1td1;
D.f = cada1f1.*CD.f;
%User Line: D = dynamic_pressure.*S.*CD;
cada1f1dV = S.*dynamic_pressure.dV;
cada1f1 = dynamic_pressure.f*S;
cada1tf1 = CL.f(:,Gator1Data.Index3);
cada1td1 = zeros(size(cada1f1dV,1),3);
cada1td1(:,Gator1Data.Index4) = cada1tf1.*cada1f1dV;
cada1td1(:,3) = cada1td1(:,3) + cada1f1.*CL.dV;
L.dV = cada1td1;
L.f = cada1f1.*CL.f;
%User Line: L = dynamic_pressure.*S.*CL;
cada1f1dV = 2.*rad.f.^(2-1).*rad.dV;
cada1f1 = rad.f.^2;
gravity.dV = -mu./cada1f1.^2.*cada1f1dV;
gravity.f = mu./cada1f1;
%User Line: gravity = mu./rad.^2;
bank.f =  0;
%User Line: bank = 0;
sgam.dV = cos(gam.f).*gam.dV;
sgam.f = sin(gam.f);
%User Line: sgam = sin(gam);
cgam.dV = -sin(gam.f).*gam.dV;
cgam.f = cos(gam.f);
%User Line: cgam = cos(gam);
cbank.f = cos(bank.f);
%User Line: cbank = cos(bank);
%User Line: % New equations needed
cada1f1 = y*R;
cada1f2 = cada1f1*T0;
va.f = sqrt(cada1f2);
%User Line: va = sqrt(y*R*T0);
M.dV = v.dV./va.f;
M.f = v.f/va.f;
%User Line: M = v./va;
cada1f1dV = -h.dV;
cada1f1 = uminus(h.f);
cada1f2dV = cada1f1dV./H;
cada1f2 = cada1f1/H;
cada1f3dV = exp(cada1f2).*cada1f2dV;
cada1f3 = exp(cada1f2);
rho.dV = rho0.*cada1f3dV;
rho.f = rho0*cada1f3;
%User Line: rho = rho0*exp(-h/H);
cada1td1 = zeros(size(rho.dV,1),2);
cada1td1(:,1) = A.f.*rho.dV;
cada1td1(:,2) = cada1td1(:,2) + rho.f.*A.dV;
cada1f1dV = cada1td1;
cada1f1 = rho.f.*A.f;
cada1tf1 = v.f(:,Gator1Data.Index5);
cada1td1 = zeros(size(cada1f1dV,1),3);
cada1td1(:,Gator1Data.Index6) = cada1tf1.*cada1f1dV;
cada1td1(:,2) = cada1td1(:,2) + cada1f1.*v.dV;
madot.dV = cada1td1;
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
cada1f2dV = 2.*M.f.^(2-1).*M.dV;
cada1f2 = M.f.^2;
cada1f3dV = cada1f1.*cada1f2dV;
cada1f3 = cada1f1*cada1f2;
cada1f4dV = cada1f3dV./2;
cada1f4 = cada1f3/2;
Tau_r.dV = cada1f4dV;
Tau_r.f = 1 + cada1f4;
%User Line: Tau_r = (1+(y-1).*M.^2/2);
%User Line: %Tau_lam = Tau_r +(Tau_lam - Tau_r).*(atan(qscale*(q-q1))-atan(qscale*(q-q2)))/pi;
cada1f1 = cp*T0;
cada1f2 = cada1f1/hpr;
cada1f3dV = -Tau_r.dV;
cada1f3 = Tau_lam.f - Tau_r.f;
f.dV = cada1f2.*cada1f3dV;
f.f = cada1f2*cada1f3;
%User Line: f = cp*T0/hpr*(Tau_lam-Tau_r);
cada1f1dV = -Tau_lam.f./Tau_r.f.^2.*Tau_r.dV;
cada1f1 = Tau_lam.f./Tau_r.f;
cada1f2dV = (1/2)./sqrt(cada1f1).*cada1f1dV;
cada1f2dV(cada1f1 == 0 & cada1f1dV == 0) = 0;
cada1f2 = sqrt(cada1f1);
cada1f3dV = cada1f2dV;
cada1f3 = cada1f2 - 1;
cada1td1 = cada1f3.*v.dV;
cada1td1 = cada1td1 + v.f.*cada1f3dV;
ST.dV = cada1td1;
ST.f = v.f.*cada1f3;
%User Line: ST = v.*(sqrt(Tau_lam./Tau_r)-1);
cada1f1dV = g0.*f.dV;
cada1f1 = f.f*g0;
cada1td1 = ST.dV./cada1f1;
cada1td1 = cada1td1 + -ST.f./cada1f1.^2.*cada1f1dV;
Isp.dV = cada1td1;
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
cada1tf1 = ST.f(:,Gator1Data.Index7);
cada1td1 = cada1tf1.*madot.dV;
cada1td1(:,2) = cada1td1(:,2) + madot.f.*ST.dV;
T.dV = cada1td1;
T.f = madot.f.*ST.f;
%User Line: T = madot.*ST;
%User Line: % T = T.*(atan(1e17*(M-M1))+atan(1e17*(M-M2)))/pi;
%User Line: % Tmax = input.auxdata.Tmax;
%User Line: % T = Tmax*A;
%User Line: % path = dynamic_pressure;
%User Line: % output.path = path;
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,1) = sgam.f.*v.dV;
cada1td1(:,2) = cada1td1(:,2) + v.f.*sgam.dV;
hdot.dV = cada1td1;
hdot.f = v.f.*sgam.f;
%User Line: hdot   = v.*sgam;
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,1) = cgam.f.*v.dV;
cada1td1(:,2) = cada1td1(:,2) + v.f.*cgam.dV;
cada1f1dV = cada1td1;
cada1f1 = v.f.*cgam.f;
cada1tf1 = rad.f(:,Gator1Data.Index8);
cada1td1 = zeros(size(cada1f1dV,1),3);
cada1td1(:,Gator1Data.Index9) = cada1f1dV./cada1tf1;
cada1td1(:,1) = cada1td1(:,1) + -cada1f1./rad.f.^2.*rad.dV;
thettadot.dV = cada1td1;
thettadot.f = cada1f1./rad.f;
%User Line: thettadot = v.*cgam./rad;
cada1f1dV = -sin(aoa.f).*aoa.dV;
cada1f1 = cos(aoa.f);
cada1tf1 = cada1f1(:,Gator1Data.Index10);
cada1td1 = zeros(size(T.dV,1),4);
cada1td1(:,Gator1Data.Index11) = cada1tf1.*T.dV;
cada1td1(:,3) = cada1td1(:,3) + T.f.*cada1f1dV;
cada1f2dV = cada1td1;
cada1f2 = T.f.*cada1f1;
cada1td1 = cada1f2dV;
cada1td1(:,Gator1Data.Index12) = cada1td1(:,Gator1Data.Index12) + -D.dV;
cada1f3dV = cada1td1;
cada1f3 = cada1f2 - D.f;
cada1tf1 = mass.f(:,Gator1Data.Index13);
cada1td1 = zeros(size(cada1f3dV,1),5);
cada1td1(:,Gator1Data.Index14) = cada1f3dV./cada1tf1;
cada1td1(:,3) = cada1td1(:,3) + -cada1f3./mass.f.^2.*mass.dV;
cada1f4dV = cada1td1;
cada1f4 = cada1f3./mass.f;
cada1td1 = zeros(size(gravity.dV,1),2);
cada1td1(:,1) = sgam.f.*gravity.dV;
cada1td1(:,2) = cada1td1(:,2) + gravity.f.*sgam.dV;
cada1f5dV = cada1td1;
cada1f5 = gravity.f.*sgam.f;
cada1td1 = zeros(size(cada1f4dV,1),6);
cada1td1(:,Gator1Data.Index15) = cada1f4dV;
cada1td1(:,Gator1Data.Index16) = cada1td1(:,Gator1Data.Index16) + -cada1f5dV;
vdot.dV = cada1td1;
vdot.f = cada1f4 - cada1f5;
%User Line: vdot = (T.*cos(aoa)-D)./mass-gravity.*sgam;
cada1f1dV = cos(aoa.f).*aoa.dV;
cada1f1 = sin(aoa.f);
cada1tf1 = cada1f1(:,Gator1Data.Index17);
cada1td1 = zeros(size(T.dV,1),4);
cada1td1(:,Gator1Data.Index18) = cada1tf1.*T.dV;
cada1td1(:,3) = cada1td1(:,3) + T.f.*cada1f1dV;
cada1f2dV = cada1td1;
cada1f2 = T.f.*cada1f1;
cada1f3dV = cbank.f.*L.dV;
cada1f3 = L.f*cbank.f;
cada1td1 = cada1f2dV;
cada1td1(:,Gator1Data.Index19) = cada1td1(:,Gator1Data.Index19) + cada1f3dV;
cada1f4dV = cada1td1;
cada1f4 = cada1f2 + cada1f3;
cada1tf1 = mass.f(:,Gator1Data.Index20);
cada1td1 = zeros(size(cada1f4dV,1),5);
cada1td1(:,Gator1Data.Index21) = cada1f4dV./cada1tf1;
cada1td1(:,3) = cada1td1(:,3) + -cada1f4./mass.f.^2.*mass.dV;
cada1f5dV = cada1td1;
cada1f5 = cada1f4./mass.f;
cada1f6dV = 2.*v.f.^(2-1).*v.dV;
cada1f6 = v.f.^2;
cada1td1 = zeros(size(cada1f6dV,1),2);
cada1td1(:,2) = cada1f6dV./rad.f;
cada1td1(:,1) = cada1td1(:,1) + -cada1f6./rad.f.^2.*rad.dV;
cada1f7dV = cada1td1;
cada1f7 = cada1f6./rad.f;
cada1td1 = zeros(size(gravity.dV,1),2);
cada1td1(:,1) = gravity.dV;
cada1td1 = cada1td1 + -cada1f7dV;
cada1f8dV = cada1td1;
cada1f8 = gravity.f - cada1f7;
cada1td1 = zeros(size(cgam.dV,1),3);
cada1td1(:,3) = cada1f8.*cgam.dV;
cada1tf1 = cgam.f(:,Gator1Data.Index22);
cada1td1(:,Gator1Data.Index23) = cada1td1(:,Gator1Data.Index23) + cada1tf1.*cada1f8dV;
cada1f9dV = cada1td1;
cada1f9 = cgam.f.*cada1f8;
cada1td1 = zeros(size(cada1f5dV,1),6);
cada1td1(:,Gator1Data.Index24) = cada1f5dV;
cada1td1(:,Gator1Data.Index25) = cada1td1(:,Gator1Data.Index25) + -cada1f9dV;
cada1f10dV = cada1td1;
cada1f10 = cada1f5 - cada1f9;
cada1tf1 = v.f(:,Gator1Data.Index26);
cada1td1 = cada1f10dV./cada1tf1;
cada1td1(:,2) = cada1td1(:,2) + -cada1f10./v.f.^2.*v.dV;
gamdot.dV = cada1td1;
gamdot.f = cada1f10./v.f;
%User Line: gamdot   = ((T.*sin(aoa)+L.*cbank)./mass-cgam.*(gravity-v.^2./rad))./v;
cada1f1dV = -madot.dV;
cada1f1 = uminus(madot.f);
cada1tf1 = f.f(:,Gator1Data.Index27);
cada1td1 = cada1tf1.*cada1f1dV;
cada1td1(:,2) = cada1td1(:,2) + cada1f1.*f.dV;
massdot.dV = cada1td1;
massdot.f = cada1f1.*f.f;
%User Line: massdot = -madot.*f;
aoadot.dV = alfadot.dV; aoadot.f = alfadot.f;
%User Line: aoadot = alfadot;
cada1td1 = zeros(size(hdot.f,1),21);
cada1td1(:,Gator1Data.Index28) = hdot.dV;
cada1td1(:,Gator1Data.Index29) = thettadot.dV;
cada1td1(:,Gator1Data.Index30) = vdot.dV;
cada1td1(:,Gator1Data.Index31) = gamdot.dV;
cada1td1(:,Gator1Data.Index32) = massdot.dV;
cada1td1(:,18) = aoadot.dV;
output.dynamics.dV = cada1td1;
output.dynamics.f = [hdot.f thettadot.f vdot.f gamdot.f massdot.f aoadot.f];
%User Line: output.dynamics  = [hdot, thettadot, vdot, gamdot, massdot, aoadot];
%User Line: %------------------------------------%
%User Line: % End File:  goddardRocketContinuous %
%User Line: %------------------------------------%
output.dynamics.dV_size = [6,9];
output.dynamics.dV_location = Gator1Data.Index33;
end


function ADiGator_LoadData()
global ADiGator_Thrust_real_updated_ContinuousADiGatorGrd
ADiGator_Thrust_real_updated_ContinuousADiGatorGrd = load('Thrust_real_updated_ContinuousADiGatorGrd.mat');
return
end