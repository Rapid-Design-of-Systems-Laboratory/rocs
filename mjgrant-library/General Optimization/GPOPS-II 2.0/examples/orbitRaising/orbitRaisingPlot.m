%------------------------------%
% Extract Solution from Output %
%------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
r         = solution.phase(1).state(:,1);
theta     = solution.phase(1).state(:,2);
vr        = solution.phase(1).state(:,3);
vtheta    = solution.phase(1).state(:,4);
u1        = solution.phase(1).control(:,1);
u2        = solution.phase(1).control(:,2);
alpha     = unwrap(atan2(u1,u2))*180/pi;

state = [r, theta, vr, vtheta];
control = [u1, u2];
%---------------%
% Plot Solution %
%---------------%
figure(1)
pp = plot(time,state,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(r(t),\theta(t),v_r(t),v_\theta(t))$','Interpreter','LaTeX');
ll = legend('$r(t)$','$\theta(t)$','$v_r(t)$','$v_\theta(t)$','Location','Northwest');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 orbitRaisingState.eps;
print -dpng orbitRaisingState.png;

figure(2)
pp = plot(time,r,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
%print -depsc2 orbitRaisingRadius.eps
print -dpng orbitRaisingRadius.png

figure(3)
pp = plot(time,theta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
%print -depsc2 orbitRaisingTheta.eps
print -dpng orbitRaisingTheta.png

figure(4)
pp = plot(time,vr,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
%print -depsc2 orbitRaisingVr.eps
print -dpng orbitRaisingVr.png

figure(5)
pp = plot(time,vtheta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on
%print -depsc2 orbitRaisingVtheta.eps
print -dpng orbitRaisingVtheta.png

figure(6)
pp = plot(time,control,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(u_1(t),u_2(t))$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 orbitRaisingControl.eps;
print -dpng orbitRaisingControl.png;

figure(7)
pp = plot(time,alpha,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\alpha(t)=\displaystyle\tan^{-1}\left(\frac{u_1(t)}{u_2(t)}\right)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 orbitRaisingAlpha.eps;
print -dpng orbitRaisingAlpha.png;
