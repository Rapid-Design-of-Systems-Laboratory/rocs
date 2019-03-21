% Extract Solution.                       
solution = output.result.solution;
time = solution.phase.time;
state = solution.phase.state;
control = solution.phase.control;
realControl = [control(:,1)-control(:,2),control(:,3)-control(:,4)];
for i=1:length(output.meshhistory);
  mesh(i).points = [0 cumsum(output.meshhistory(i).result.setup.mesh.phase.fraction)];
  mesh(i).iteration = i*ones(size(mesh(i).points));
end;

% Plot Solution.
figure(1);
pp = plot(time,state(:,1),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$x(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng Time-vs-xr.png

figure(2);
pp = plot(time,state(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$y(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
%axis square
grid on
print -dpng Time-vs-yr.png

figure(3);
pp = plot(time,state(:,3),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\theta(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
%axis square
grid on
print -dpng Time-vs-thetar.png

figure(4);
pp = plot(time,state(:,4),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_x(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
%axis square
grid on
print -dpng Time-vs-Vxr.png

figure(5);
pp = plot(time,state(:,5),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_y(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
%axis square
grid on
print -dpng Time-vs-Vyr.png

figure(6);
pp = plot(time,state(:,6),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\omega(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
%axis square
grid on
print -dpng Time-vs-angVel.png

figure(7);
pp = plot(time,realControl(:,1),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$T_1(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
%axis square
grid on
print -dpng Time-vs-Control-1.png

figure(8);
pp = plot(time,realControl(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$T_2(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
%axis square
grid on
print -dpng Time-vs-Control-2.png

