%------------------------------%
% Extract Solution from Output %
%------------------------------%
solution = output.result.solution;

time = solution.phase.time;
down  = solution.phase.state(:,1)/1000;
alt = solution.phase.state(:,2)/1000;
v = solution.phase.state(:,3);
gam  = solution.phase.state(:,4)*180/pi;
aoa = solution.phase.control(:,1)*180/pi;
T = solution.phase.control(:,2);

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

figure(4)
plot(time,Noise,'b', 'markersize', 7, 'linewidth', 2);
yl = xlabel('Time (s)');
xl = ylabel('Noise (dB)');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16);
title('Noise History Plot','FontSize',18)
grid on

return