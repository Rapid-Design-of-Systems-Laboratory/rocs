%------------------------------%
% Extract Solution from Output %
%------------------------------%
solution = output.result.solution;
% Time
time = solution.phase(1).time;
% Phase
q1         = solution.phase(1).state(:,1);
q2         = solution.phase(1).state(:,2);
q3         = solution.phase(1).state(:,3);
w1         = solution.phase(1).state(:,4);
w2         = solution.phase(1).state(:,5);
w3         = solution.phase(1).state(:,6);
% Control
q4        = solution.phase(1).control(:,1);
u1        = solution.phase(1).control(:,2);
u2        = solution.phase(1).control(:,3);
u3        = solution.phase(1).control(:,4);

%---------------%
% Plot Solution %
%---------------%

% Plot the quaternions. ------------------------------------------------------

figure(1)
pp = plot(time,q1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$q_1(t)$','Interpreter','LaTeX');
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 q1.eps
print -dpng q1.png

figure(2)
pp = plot(time,q2,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$q_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 q2.eps
print -dpng q2.png

figure(3)
pp = plot(time,q3,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$q_3(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 q3.eps
print -dpng q3.png

figure(4)
pp = plot(time,q4,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$q_4(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on
print -depsc2 q4.eps
print -dpng q4.png

% Plot the ang vel. ------------------------------------------------------

figure(5)
pp = plot(time,w1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\omega_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 w1.eps
print -dpng w1.png

figure(6)
pp = plot(time,w2,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\omega_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 w2.eps
print -dpng w2.png

figure(7)
pp = plot(time,w3,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\omega_3(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 w3.eps
print -dpng w3.png

% Plot the Controls ------------------------------------------------------
figure(8)
pp = plot(time,u1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 u1.eps
print -dpng u1.png

figure(9)
pp = plot(time,u2,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 u2.eps
print -dpng u2.png

figure(10)
pp = plot(time,u3,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_3(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 u3.eps
print -dpng u3.png