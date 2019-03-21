solution = output.result.solution;
time  = solution.phase.time;
state = solution.phase.state;
control = solution.phase.control;

ifig = 1;
figure(ifig);
pp = plot(time,state(:,1)*1e4,'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\omega_1(t)$~(rad$\cdot$s${}^{-1}\times{10^{-4}}$)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState1.eps
print -dpng spaceStationState1.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,2)*1e4,'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\omega_2(t)$~(rad$\cdot$s${}^{-1}\times{10^{-4}}$)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState2.eps
print -dpng spaceStationState2.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,3)*1e4,'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\omega_3(t)$~(rad$\cdot$s${}^{-1}\times{10^{-4}}$)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState3.eps
print -dpng spaceStationState3.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,4),'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$r_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState4.eps
print -dpng spaceStationState4.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,5),'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$r_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState5.eps
print -dpng spaceStationState5.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,6),'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$r_3(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState6.eps
print -dpng spaceStationState6.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,7)/1000,'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$h_1(t)$~(N$\cdot$m$\cdot$s)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState7.eps
print -dpng spaceStationState7.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,8)/1000,'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$h_2(t)$~(N$\cdot$m$\cdot$s)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState8.eps
print -dpng spaceStationState8.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,state(:,9)/1000,'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$h_3(t)$~(N$\cdot$m$\cdot$s)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState9.eps
print -dpng spaceStationState9.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,control(:,1),'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$~(N$\cdot$m)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationControl1.eps
print -dpng spaceStationControl1.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,control(:,2),'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$u_2(t)$~(N$\cdot$m)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationControl2.eps
print -dpng spaceStationControl2.png

ifig = ifig+1;
figure(ifig);
pp = plot(time,control(:,3),'-o');
set(pp,'LineWidth',1.25,'MarkerSize',6);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$u_3(t)$~(N$\cdot$m)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationControl3.eps
print -dpng spaceStationControl3.png
