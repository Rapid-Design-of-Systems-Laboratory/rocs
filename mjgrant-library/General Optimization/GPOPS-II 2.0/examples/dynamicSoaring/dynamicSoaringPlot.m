t = solution.phase(1).time;
x = solution.phase(1).state(:,1);
y = solution.phase(1).state(:,2);
z = solution.phase(1).state(:,3);
v = solution.phase(1).state(:,4);
gamma = solution.phase(1).state(:,5);
psi = solution.phase(1).state(:,6);
CL = solution.phase(1).control(:,1);
bank = solution.phase(1).control(:,2);

figure(1);
pp = plot3(x,y,z,'-o');
xl = xlabel('$x(t)$~(ft)','Interpreter','LaTeX');
yl = ylabel('$y(t)$~(ft)','Interpreter','LaTeX');
zl = zlabel('$z(t)$~(ft)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(zl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -dpng dynamicSoaringXYZ.png

figure(2);
pp = plot(t,v,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$v(t)$~(ft$\cdot$s${}^{-1}$)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'YTick',40:40:240,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -dpng dynamicSoaringV.png

figure(3);
pp = plot(t,gamma*180/pi,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\gamma(t)$~(deg)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -dpng dynamicSoaringGamma.png

figure(4);
pp = plot(t,psi*180/pi,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\psi(t)$~(deg)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -dpng dynamicSoaringPsi.png

figure(5);
pp = plot(t,CL,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$C_L(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -dpng dynamicSoaringCL.png

figure(6);
pp = plot(t,bank*180/pi,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$\sigma(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -dpng dynamicSoaringSigma.png
