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
pp = plot(time,state,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$x(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng hyperSensitiveState.png

figure(2);
pp = plot(time,control,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng hyperSensitiveControl.png

figure(3);
for i=1:length(mesh);
  pp = plot(mesh(i).meshPoints*tf,mesh(i).iteration,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Mesh Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16,'FontName','Times');
axis([0 10000 1 7]);
grid on;
print -dpng hyperSensitiveMeshHistory.png

figure(4);
for i=1:length(mesh);
  pp = plot(mesh(i).time,mesh(i).iterationTime,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Collocation Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16,'FontName','Times');
axis([0 10000 1 7]);
grid on;
print -dpng hyperSensitiveCollocHistory.png
