%-------------------------------------------------------------------%
%-------------------------- Plot Solution --------------------------%
%-------------------------------------------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
S    = solution.phase(1).state(:,1);
L1   = solution.phase(1).state(:,2);
L2   = solution.phase(1).state(:,3);
I1   = solution.phase(1).state(:,4);
I2   = solution.phase(1).state(:,5);
T    = solution.phase(1).state(:,6);
u1   = solution.phase(1).control(:,1);
u2   = solution.phase(1).control(:,2);

figure(1);
pp = plot(time,S/1000,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$S(t)\times 10^3$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisS.eps
print -dpng tuberculosisS.png

figure(2);
pp = plot(time,L1/1000,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$L_1(t)\times 10^3$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisL1.eps
print -dpng tuberculosisL1.png

figure(3);
pp = plot(time,L2/1000,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$L_2(t)\times 10^3$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisL2.eps
print -dpng tuberculosisL2.png

figure(4);
pp = plot(time,I1,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$I_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisI1.eps
print -dpng tuberculosisI1.png

figure(5);
pp = plot(time,I2/1000,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$I_2(t)\times 10^3$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisI2.eps
print -dpng tuberculosisI2.png

figure(6);
pp = plot(time,T,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$T(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisT.eps
print -dpng tuberculosisT.png

figure(7);
pp = plot(time,u1,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$U_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisU1.eps
print -dpng tuberculosisU1.png

figure(8);
pp = plot(time,u2,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$U_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on;
print -depsc2 tuberculosisU2.eps
print -dpng tuberculosisU2.png
