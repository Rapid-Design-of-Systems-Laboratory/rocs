%------------------------------%
% Extract Solution from Output %
%------------------------------%

solution = output.result.solution;
save solution
close all; clc
% 
time = solution.phase.time;
h  = solution.phase.state(:,1);
lamH = solution.phase.costate(:,1);
lamTHETTA = solution.phase.costate(:,2);
lamV = solution.phase.costate(:,3);
lamGAM = solution.phase.costate(:,4);
lamMASS = solution.phase.costate(:,5);
lamAOA = solution.phase.costate(:,6);
downrange = 6378*solution.phase.state(:,2);
v = solution.phase.state(:,3);
gam  = solution.phase.state(:,4);
mass = solution.phase.state(:,5);
aoa = solution.phase.state(:,6)*180/pi;
aoadot = solution.phase.control(:,1)*180/pi;
A = solution.phase.control(:,2);
fixed = solution;
save fixed
% Earth
Re   = 6378000;                     % Equatorial Radius of Earth (m)
H    = 7500;                        % Density Scale Height (m)
rho0 = 1.2;              % Sea Level Atmospheric Density (kg/m^3)
mu   = 3.986e5*1e9;           % Earth Gravitational Parameter (m^3/s^2)
g0 = 9.81; 

% Isp = 1600;
Tmax = 1800000;

% New parameters

y  = 1.4;
R  = 287.058;
T0  = 230;
T4  = 1600;
cp  = 1004;
hpr  = 43903250;
Mc  = 3;
M1  = 4;
M2  = 9.75;
q1  = 40000;
q2  = 150000;
T1  = 0;
rho = rho0*exp(-h/H); % Exponential Atmospheric Density [kg/m^3]
q = 0.5*rho.*v.^2; % Dynamic Pressure
% New equations needed
va = sqrt(y*R*T0);
M = v./va;

madot = rho.*A.*v; % air flow rate
Tau_lam = T4*(1+(y-1)*Mc^2/2)/T0;
Tau_r = (1+(y-1).*M.^2/2);
% Tau_lam = Tau_r +(Tau_lam - Tau_r).*(atan(1e3*(dynamic_pressure-q1))-atan(1e3*(dynamic_pressure-q2)))/pi;
f = cp*T0/hpr*(Tau_lam-Tau_r); % fuel ratio
ST = v.*(sqrt(Tau_lam./Tau_r)-1);
Isp = ST./(f*g0);
% for i = 1:1:length(ST)
%     % M_info = M(i);
%     q_info = q(i);
%     if (q_info<50000 || q_info>150000) % M_info<4 || M_info>9.75 || 
%         T(i) = 0;
%     else
%       T(i) = madot(i)*ST(i);  
%     end
% end
T = madot.*ST;
aoa_rad = aoa*pi/180;
aoadot_rad = aoadot*pi/180;
cl   = [0.1758 10.305];
cd   = [0.26943 -0.4113 18.231];
cd0      = cd(1);
cd1      = cd(2);
cd2      = cd(3);
cl0      = cl(1);
cl1      = cl(2);
CD       = cd0+cd1*aoa_rad+cd2*aoa_rad.^2;
CL       = cl0+cl1*aoa_rad;
S = 0.35;
D = q.*S.*CD;
L = q.*S.*CL;

sgam = sin(gam);
rad = h + Re;
gravity = mu./rad.^2;
cgam = cos(gam);
hdot   = v.*sgam;
thettadot = v.*cgam./rad;
vdot = (T.*cos(aoa_rad)-D)./mass-gravity.*sgam;
gamdot   = ((T.*sin(aoa_rad)+L)./mass-cgam.*(gravity-v.^2./rad))./v;
massdot = -madot.*f;
alfadot = aoadot_rad;

Ham = lamH.*hdot + lamTHETTA.*thettadot + lamV.*vdot + lamGAM.*gamdot + lamMASS.*massdot + lamAOA.*alfadot;
gam = gam.*180./pi;
h = h./1000;
%----------%
% Plot Solution %
%---------------%
figure(1)
subplot(1,2,1)
plot(v/1000,h,'b', 'markersize', 7, 'linewidth', 2);
xl = xlabel('Velocity (km/s)');
yl = ylabel('Altitude (km)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Energy Plot','FontSize',18)
grid on

subplot(1,2,2)
plot(time,gam,'b', 'markersize', 7, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('Flight Path Angle (deg)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Flight Path Angle History Plot','FontSize',18)
grid on

figure(2)
subplot(1,2,1)
plot(downrange,h,'k', 'markersize', 7, 'linewidth', 2);
xl = xlabel('Downrange (km)');
yl = ylabel('Altitude (km)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
grid on
title('Trajectory Plot','FontSize',18)

subplot(1,2,2)
plot(time,mass,'k', 'markersize', 7, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('Mass (kg)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
grid on
title('Mass time history','FontSize',18)



figure(3)
subplot(1,3,1)
plot(time,aoa,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Angle of Attack (deg)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('AoA History Plot','FontSize',18)
grid on

subplot(1,3,2)
plot(time,A,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Area (m^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Thrust Control History Plot','FontSize',18)
grid on

subplot(1,3,3)
plot(time,aoadot,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Rate of AoA (deg/s)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Rate of AoA Control History Plot','FontSize',18)
grid on


figure(4)
subplot(1,3,1)
plot(time,q,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Dynamic Pressure (N/m^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Dynamic Pressure History Plot','FontSize',18)
grid on

subplot(1,3,2)
plot(time,M,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('M ()');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Mach History Plot','FontSize',18)
grid on

subplot(1,3,3)
plot(time,Isp,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Isp (s)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Isp History Plot','FontSize',18)
grid on

figure(5)
subplot(1,3,1)
plot(time,L,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Lift (N)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Lift History Plot','FontSize',18)
grid on

subplot(1,3,2)
plot(time,D,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Drag (N)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Drag History Plot','FontSize',18)
grid on

subplot(1,3,3)
plot(time,T,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Thrust (N)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Thrust History Plot','FontSize',18)
grid on

figure(6)
subplot(1,3,1)
plot(time,h,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Altitude (km)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Altitude History Plot','FontSize',18)
grid on

subplot(1,3,2)
plot(time,downrange,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Downrange (km)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Downrange History Plot','FontSize',18)
grid on
subplot(1,3,3)
plot(time,v/1000,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Velocity (m/s)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on

figure(7)
% subplot(1,2,1)
plot(time,vdot,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Vdot (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
% subplot(1,2,2)
% plot(time,(T.*cos(aoa_rad)-D)./mass,'b', 'markersize', 7, 'linewidth', 2);
% yl = xlabel('Time (s)');
% xl = ylabel('Axial Forces (N)');
% set(xl,'FontSize',16);
% set(yl,'FontSize',16);
% set(gca,'FontSize',16);
% title('Axial Force History Plot','FontSize',18)
% grid on
figure(8)
subplot(1,3,1)
plot(time,lamH,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('lamH (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
subplot(1,3,2)
plot(time,lamTHETTA,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('lamTHETTA (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
subplot(1,3,3)
plot(time,lamV,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('lamV (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
figure(9)
subplot(1,3,1)
plot(time,lamGAM,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('lamGAM (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
subplot(1,3,2)
plot(time,lamMASS,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('lamMASS (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
subplot(1,3,3)
plot(time,lamAOA,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('lamAOA (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
figure(10)
plot(time,Ham,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Hamiltonian (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Hamiltonian History Plot','FontSize',18)
grid on
figure(11)
subplot(1,2,1)
plot(time,Tau_lam - Tau_r,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Tau_lam - taur (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
subplot(1,2,2)
plot(time,Tau_lam./Tau_r,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('taulam/taur (m/s^2)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Velocity History Plot','FontSize',18)
grid on
return