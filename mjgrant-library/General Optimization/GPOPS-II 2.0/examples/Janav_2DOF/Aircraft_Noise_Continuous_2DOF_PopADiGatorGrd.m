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

function output = Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd(input)
global ADiGator_Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd
if isempty(ADiGator_Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd.Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %----------------------------------------%
%User Line: % Begin File:  Aircraft Noise Continuous %
%User Line: %----------------------------------------%
%User Line: % Constants
mass = input.auxdata.mass;
%User Line: mass = input.auxdata.mass;
g0 = input.auxdata.g0;
%User Line: g0 = input.auxdata.g0;
C3 = input.auxdata.C3;
%User Line: C3 = input.auxdata.C3;
C1 = input.auxdata.C1;
%User Line: C1 = input.auxdata.C1;
C2 = input.auxdata.C2;
%User Line: C2 = input.auxdata.C2;
%User Line: % States and time
t.dV = input.phase.time.dV; t.f = input.phase.time.f;
%User Line: t = input.phase.time;
state.dV = input.phase.state.dV; state.f = input.phase.state.f;
%User Line: state = input.phase.state;
x.dV = state.dV(:,1);
x.f = state.f(:,1);
%User Line: x = state(:,1);
z.dV = state.dV(:,2);
z.f = state.f(:,2);
%User Line: z = state(:,2);
v.dV = state.dV(:,3);
v.f = state.f(:,3);
%User Line: v = state(:,3);
gam.dV = state.dV(:,4);
gam.f = state.f(:,4);
%User Line: gam = state(:,4);
%User Line: % Controls
aoa.dV = input.phase.control.dV(:,1);
aoa.f = input.phase.control.f(:,1);
%User Line: aoa  = input.phase.control(:,1);
T.dV = input.phase.control.dV(:,2);
T.f = input.phase.control.f(:,2);
%User Line: T = input.phase.control(:,2);
%User Line: % Important parameters
%User Line: % rho = rho0*exp(-h/H);
%User Line: % CD       = cd0+cd1*aoa+cd2*aoa.^2;
%User Line: % CL       = cl0+cl1*aoa;
%User Line: % dynamic_pressure = 0.5*rho.*v.^2;
%User Line: % D = dynamic_pressure.*S.*CD;
%User Line: % L = dynamic_pressure.*S.*CL;
L.f = mass*g0;
%User Line: L = mass.*g0;
cada1f1dV = 2.*v.f.^(2-1).*v.dV;
cada1f1 = v.f.^2;
cada1f2dV = 0.226.*cada1f1dV;
cada1f2 = 0.226*cada1f1;
cada1f3dV = 2.*v.f.^(2-1).*v.dV;
cada1f3 = v.f.^2;
cada1f4dV = -5200000./cada1f3.^2.*cada1f3dV;
cada1f4 = 5200000./cada1f3;
cada1td1 = cada1f2dV;
cada1td1 = cada1td1 + cada1f4dV;
D.dV = cada1td1;
D.f = cada1f2 + cada1f4;
%User Line: D = 0.226.*v.^2 + 5.2e6./(v.^2);
gravity = g0;
%User Line: gravity = g0;
sgam.dV = cos(gam.f).*gam.dV;
sgam.f = sin(gam.f);
%User Line: sgam = sin(gam);
cgam.dV = -sin(gam.f).*gam.dV;
cgam.f = cos(gam.f);
%User Line: cgam = cos(gam);
%User Line: % Approximate population
%User Line: % keyboard
%User Line: % for count = 1:1:length(x)
%User Line: % pop(count) = population(x(count)/1000,y(count)/1000);
%User Line: % end
%User Line: % pop = -0.093425*x.^3 + 0.017248*x.^2.*y +...
%User Line: %     0.040138*x.^2 - 0.0035582*x.*y.^2 + 0.076212*x.*y +...
%User Line: %     1.3537*x - 0.0069758*x.^3 - 0.059675*y.^2 + 0.062387*y + 2.481;
%User Line: % Below code works
%User Line: % pop = 1;
%User Line: %
%User Line: % integrand  = log(C3)+5.2*log(T)+log(cos(gam))-(log(v)+2.5*log(h+50));
%User Line: % x1 = x/1000;
%User Line: % y1 = y/1000;
%User Line: % pop = 1;
cada1f1dV = -x.dV;
cada1f1 = 15000 - x.f;
pop.dV = 2.*cada1f1dV;
pop.f = 2*cada1f1;
%User Line: pop = 2*(15000-x);
%User Line: % pop = -0.093425*x1.^3 + 0.017248.*x1.^2.*y1 + 0.74082.*x1.^2 - 0.0035582.*x1.*y1.^2 - 0.010029.*x1.*y1 - 0.59874.*x1 - 0.0069758.*y1.^3 - 0.05078.*y1.^2 - 0.020341.*y1 + 0.8075;
cada1f1dV = C3.*pop.dV;
cada1f1 = pop.f*C3;
cada1f2dV = 5.2.*T.f.^(5.2-1).*T.dV;
cada1f2 = T.f.^5.2;
cada1td1 = zeros(size(cada1f1dV,1),2);
cada1td1(:,1) = cada1f2.*cada1f1dV;
cada1td1(:,2) = cada1td1(:,2) + cada1f1.*cada1f2dV;
cada1f3dV = cada1td1;
cada1f3 = cada1f1.*cada1f2;
cada1f4dV = -sin(gam.f).*gam.dV;
cada1f4 = cos(gam.f);
cada1tf1 = cada1f4(:,Gator1Data.Index1);
cada1td1 = zeros(size(cada1f3dV,1),3);
cada1td1(:,Gator1Data.Index2) = cada1tf1.*cada1f3dV;
cada1td1(:,2) = cada1td1(:,2) + cada1f3.*cada1f4dV;
cada1f5dV = cada1td1;
cada1f5 = cada1f3.*cada1f4;
cada1f6dV = z.dV;
cada1f6 = z.f + 50;
cada1f7dV = 2.5.*cada1f6.^(2.5-1).*cada1f6dV;
cada1f7 = cada1f6.^2.5;
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,2) = cada1f7.*v.dV;
cada1td1(:,1) = cada1td1(:,1) + v.f.*cada1f7dV;
cada1f8dV = cada1td1;
cada1f8 = v.f.*cada1f7;
cada1tf1 = cada1f8(:,Gator1Data.Index3);
cada1td1 = zeros(size(cada1f5dV,1),5);
cada1td1(:,Gator1Data.Index4) = cada1f5dV./cada1tf1;
cada1tf1 = cada1f5(:,Gator1Data.Index5);
cada1tf2 = cada1f8(:,Gator1Data.Index6);
cada1td1(:,Gator1Data.Index7) = cada1td1(:,Gator1Data.Index7) + -cada1tf1./cada1tf2.^2.*cada1f8dV;
integrand.dV = cada1td1;
integrand.f = cada1f5./cada1f8;
%User Line: integrand = pop.*C3.*T.^(5.2).*cos(gam)./(v.*(z+50).^(2.5));
%User Line: % integrand = C3*T.^2./(v.*(h+50));
%User Line: % Equations of motion
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,1) = cgam.f.*v.dV;
cada1td1(:,2) = cada1td1(:,2) + v.f.*cgam.dV;
xdot.dV = cada1td1;
xdot.f = v.f.*cgam.f;
%User Line: xdot = v.*cgam;
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,1) = sgam.f.*v.dV;
cada1td1(:,2) = cada1td1(:,2) + v.f.*sgam.dV;
zdot.dV = cada1td1;
zdot.f = v.f.*sgam.f;
%User Line: zdot   = v.*sgam;
cada1f1dV = -sin(aoa.f).*aoa.dV;
cada1f1 = cos(aoa.f);
cada1td1 = zeros(size(T.dV,1),2);
cada1td1(:,2) = cada1f1.*T.dV;
cada1td1(:,1) = cada1td1(:,1) + T.f.*cada1f1dV;
cada1f2dV = cada1td1;
cada1f2 = T.f.*cada1f1;
cada1td1 = zeros(size(cada1f2dV,1),3);
cada1td1(:,Gator1Data.Index8) = cada1f2dV;
cada1td1(:,1) = cada1td1(:,1) + -D.dV;
cada1f3dV = cada1td1;
cada1f3 = cada1f2 - D.f;
cada1f4dV = cada1f3dV./mass;
cada1f4 = cada1f3/mass;
cada1f5dV = gravity.*sgam.dV;
cada1f5 = gravity*sgam.f;
cada1td1 = zeros(size(cada1f4dV,1),4);
cada1td1(:,Gator1Data.Index9) = cada1f4dV;
cada1td1(:,2) = cada1td1(:,2) + -cada1f5dV;
vdot.dV = cada1td1;
vdot.f = cada1f4 - cada1f5;
%User Line: vdot = (T.*cos(aoa)-D)./mass-gravity.*sgam;
cada1f1dV = cos(aoa.f).*aoa.dV;
cada1f1 = sin(aoa.f);
cada1td1 = zeros(size(T.dV,1),2);
cada1td1(:,2) = cada1f1.*T.dV;
cada1td1(:,1) = cada1td1(:,1) + T.f.*cada1f1dV;
cada1f2dV = cada1td1;
cada1f2 = T.f.*cada1f1;
cada1f3dV = cada1f2dV;
cada1f3 = cada1f2 + L.f;
cada1f4dV = cada1f3dV./mass;
cada1f4 = cada1f3/mass;
cada1f5dV = gravity.*cgam.dV;
cada1f5 = cgam.f*gravity;
cada1td1 = zeros(size(cada1f4dV,1),3);
cada1td1(:,Gator1Data.Index10) = cada1f4dV;
cada1td1(:,1) = cada1td1(:,1) + -cada1f5dV;
cada1f6dV = cada1td1;
cada1f6 = cada1f4 - cada1f5;
cada1tf1 = v.f(:,Gator1Data.Index11);
cada1td1 = zeros(size(cada1f6dV,1),4);
cada1td1(:,Gator1Data.Index12) = cada1f6dV./cada1tf1;
cada1td1(:,1) = cada1td1(:,1) + -cada1f6./v.f.^2.*v.dV;
gamdot.dV = cada1td1;
gamdot.f = cada1f6./v.f;
%User Line: gamdot   = ((T.*sin(aoa)+L)./mass-cgam.*gravity)./v;
cada1td1 = zeros(size(xdot.f,1),12);
cada1td1(:,Gator1Data.Index13) = xdot.dV;
cada1td1(:,Gator1Data.Index14) = zdot.dV;
cada1td1(:,Gator1Data.Index15) = vdot.dV;
cada1td1(:,Gator1Data.Index16) = gamdot.dV;
output.dynamics.dV = cada1td1;
output.dynamics.f = [xdot.f zdot.f vdot.f gamdot.f];
%User Line: output.dynamics  = [xdot, zdot ,vdot, gamdot];
output.integrand.dV = integrand.dV; output.integrand.f = integrand.f;
%User Line: output.integrand = integrand;
%User Line: %------------------------------------%
%User Line: % End File:  goddardRocketContinuous %
%User Line: %------------------------------------%
output.dynamics.dV_size = [4,7];
output.dynamics.dV_location = Gator1Data.Index17;
output.integrand.dV_size = 7;
output.integrand.dV_location = Gator1Data.Index18;
end


function ADiGator_LoadData()
global ADiGator_Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd
ADiGator_Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd = load('Aircraft_Noise_Continuous_2DOF_PopADiGatorGrd.mat');
return
end