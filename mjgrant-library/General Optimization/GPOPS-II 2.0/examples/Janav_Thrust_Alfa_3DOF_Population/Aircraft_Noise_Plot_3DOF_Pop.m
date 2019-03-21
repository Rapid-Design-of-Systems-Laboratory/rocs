%------------------------------%
% Extract Solution from Output %
%------------------------------%
solution = output.result.solution;
close all; clc
% 
time = solution.phase.time;
down  = solution.phase.state(:,1)/1000;
cross = solution.phase.state(:,2)/1000;
alt = solution.phase.state(:,3)/1000;
v = solution.phase.state(:,4);
psii = solution.phase.state(:,5)*180/pi;
gam  = solution.phase.state(:,6)*180/pi;
bank = solution.phase.control(:,1)*180/pi;
aoa = solution.phase.control(:,2)*180/pi;
T = solution.phase.control(:,3);
% save('Janav_Gam_Max_0_Min_minus8.mat.mat')
% save('Janav_Gam_Up.mat')

% Code for Loading Data File
% load('Janav_Gam_Up.mat')
% load('Janav_Gam_Max_0_Min_minus8.mat')

Noise = 10*log10(((18.73*T.^(5.2).*cos(gam.*pi./180))./(v.*(alt.*1000+50).^(2.5))));
save('results_gpops.mat')
%---------------%
% Plot Solution %
%---------------%
figure(1)
subplot(1,2,1)
plot(v/1000,alt,'b', 'markersize', 7, 'linewidth', 2);
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
plot(down,alt,'k', 'markersize', 7, 'linewidth', 2);
xl = xlabel('Downrange (km)');
yl = ylabel('Altitude (km)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
grid on
title('Trajectory Plot','FontSize',18)

figure(3)
plot(cross,alt,'k', 'markersize', 7, 'linewidth', 2);
xl = xlabel('Crossrange (km)');
yl = ylabel('Altitude (km)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
grid on
title('Trajectory Plot','FontSize',18)

figure(4)
subplot(1,3,1)
plot(time,bank,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Bank Angle (deg)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Bank Angle History Plot','FontSize',18)
grid on

subplot(1,3,2)
plot(time,aoa,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Angle of Attack (deg)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('AoA History Plot','FontSize',18)
grid on

subplot(1,3,3)
plot(time,T,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Thrust (N)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Thrust Control History Plot','FontSize',18)
grid on

figure(5)
plot(time,Noise,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Noise (dB)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Noise History Plot','FontSize',18)
grid on

return