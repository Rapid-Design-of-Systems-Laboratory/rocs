%-------------------------------------------------------------------------%
%                             Extract Solution                            %
%-------------------------------------------------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
state = solution.phase(1).state;
control = solution.phase(1).control;
for i=1:length(output.meshhistory);
  mesh(i).meshPoints = [0 cumsum(output.meshhistory(i).result.setup.mesh.phase.fraction)];
  mesh(i).time =  output.meshhistory(i).result.solution.phase.time;
  mesh(i).iteration = i*ones(size(mesh(i).meshPoints));
  mesh(i).iterationTime = i*ones(size(mesh(i).time));
end;

%-------------------------------------------------------------------------%
%                              Plot Solution                              %
%-------------------------------------------------------------------------%
figure(1);
pp = plot(time,state(:,1),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$y_1(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-Y1.eps
print -dpng robotArm-Y1.png

figure(2);
pp = plot(time,state(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$y_2(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-Y2.eps
print -dpng robotArm-Y2.png

figure(3);
pp = plot(time,state(:,3),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$y_3(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-Y3.eps
print -dpng robotArm-Y3.png

figure(4);
pp = plot(time,state(:,4),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$y_4(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-Y4.eps
print -dpng robotArm-Y4.png

figure(5);
pp = plot(time,state(:,5),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$y_5(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-Y5.eps
print -dpng robotArm-Y5.png

figure(6);
pp = plot(time,state(:,6),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$y_6(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-Y6.eps
print -dpng robotArm-Y6.png

figure(7);
pp = plot(time,control(:,1),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-U1.eps
print -dpng robotArm-U1.png

figure(8);
pp = plot(time,control(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_2(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-U2.eps
print -dpng robotArm-U2.png

figure(9);
pp = plot(time,control(:,3),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_3(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 robotArm-U3.eps
print -dpng robotArm-U3.png

figure(10);
tf = time(end);
for i=1:length(mesh);
  pp = plot(mesh(i).meshPoints*tf,mesh(i).iteration,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Mesh Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16);
axis([0 tf 1 mesh(end).iteration(1)]);
grid on;
print -depsc2 robotArmMeshRefinement.eps
print -dpng robotArmMeshRefinement.png

figure(5);
for i=1:length(mesh);
  pp = plot(mesh(i).time,mesh(i).iterationTime,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Collocation Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16);
axis([0 tf 1 mesh(end).iteration(1)]);
grid on;
print -depsc2 robotArmMeshRefinementTime.eps
print -dpng robotArmMeshRefinementTime.png


figure
% Location of Discontinuities
discontinuities = [2.286 2.797 4.557 6.343 6.855];

marks = {'kx','bo','gd','rs','m^','y*','cv'};
% Interval Divisions
divisions{1} = [];
divisions{2} = [2.1542 2.3065 2.7636 4.6029 6.3479 6.3901 6.8771 6.9995];
divisions{3} = [2.2811 2.3022 2.7834 4.5739 6.848 6.872];
divisions{4} = [2.286 2.8567 6.3455 6.8536];
divisions{5} = [2.7885 6.3458 6.3477];
divisions{6} = [2.7933];
divisions{7} = [2.8109];
divisions{8} = [2.7945];
divisions{9} = [2.7956];
pp= plot(discontinuities,ones(size(discontinuities)),marks{1},divisions{2},ones(size(divisions{2})),marks{2},divisions{3},ones(size(divisions{3})),marks{3},divisions{4},ones(size(divisions{4})),marks{4},divisions{5},ones(size(divisions{5})),marks{5},divisions{6},ones(size(divisions{6})),marks{6},divisions{7},ones(size(divisions{7})),marks{7});
legend('Discontinuity Location','Mesh 2','Mesh 3','Mesh 4','Mesh 5','Mesh 6','Mesh 7');
