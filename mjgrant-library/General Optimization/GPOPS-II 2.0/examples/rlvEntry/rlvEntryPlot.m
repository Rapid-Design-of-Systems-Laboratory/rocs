%------------------------------%
% Extract Solution from Output %
%------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
altitude  = (solution.phase(1).state(:,1)-auxdata.Re)/1000;
longitude = solution.phase(1).state(:,2)*180/pi;
latitude  = solution.phase(1).state(:,3)*180/pi;
speed     = solution.phase(1).state(:,4)/1000;
fpa       = solution.phase(1).state(:,5)*180/pi;
azimuth   = solution.phase(1).state(:,6)*180/pi;
aoa       = solution.phase(1).control(:,1)*180/pi;
bank      = solution.phase(1).control(:,2)*180/pi;

%---------------%
% Plot Solution %
%---------------%
figure(1)
pp = plot(time,altitude,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$h(t)$~(km)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 rlvEntryAltitude.eps
print -dpng rlvEntryAltitude.png

figure(2)
plot(longitude,latitude,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$\theta(t)$~(deg)','Interpreter','LaTeX');
yl = ylabel('$\phi(t)$~(deg)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 rlvEntryLonLat.eps
print -dpng rlvEntryLonLat.png

figure(3)
plot(time,speed,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$v(t)$~(km$\cdot$s${}^{-1}$)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 rlvEntrySpeed.eps
print -dpng rlvEntrySpeed.png

figure(4)
plot(time,fpa,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\gamma(t)$~(deg)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 rlvEntryFlightPathAngle.eps
print -dpng rlvEntryFlightPathAngle.png

figure(5)
plot(time,azimuth,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\psi(t)$~(deg)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 rlvEntryAzimuthAngle.eps
print -dpng rlvEntryAzimuthAngle.png

figure(6)
plot(time,aoa,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\alpha(t)$~(deg)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16,'YTick',[16.5 17 17.5]);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 rlvEntryAngleofAttack.eps
print -dpng rlvEntryAngleofAttack.png

figure(7)
plot(time,bank,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\sigma(t)$~(deg)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16,'YTick',[16.5 17 17.5]);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 rlvEntryBankAngle.eps
print -dpng rlvEntryBankAngle.png

for i=1:length(output.meshhistory),
  mesh(i).cumfraction = [0 cumsum(output.meshhistory(1).result.setup.mesh.phase(1).fraction)];
  mesh(i).sizecumfraction = i*ones(size(mesh(i).cumfraction));
end;
marks = {'bd','gs','r^','cv','mo','kh'};
figure(8);
pp = plot(mesh(1).cumfraction,mesh(1).sizecumfraction,marks{1},mesh(2).cumfraction,mesh(2).sizecumfraction,marks{2},mesh(3).cumfraction,mesh(3).sizecumfraction,marks{3});
yl = xlabel('Location of Mesh Points (Percent)');
xl = ylabel('Mesh Refinement Iteration Number');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(gca,'YTick',0:1:length(output.meshhistory));
axis([0 1 0 3]);
grid on
print -depsc2 rlvEntryMeshHistory.png
print -dpng rlvEntryMeshHistory.png

