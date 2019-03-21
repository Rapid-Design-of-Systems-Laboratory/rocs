%-------------------------------------------------------------------%
%--------------------------- Plot Solution -------------------------%
%-------------------------------------------------------------------%
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

tf = solution.phase.time(end);

figure(1);
pp = plot(solution.phase(1).time,solution.phase(1).state(:,1),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$p(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on;
print -depsc2 tumorAntiAngiogenesisState1.eps
print -dpng tumorAntiAngiogenesisState1.png

figure(2);
pp = plot(solution.phase(1).time,solution.phase(1).state(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$q(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on;
print -depsc2 tumorAntiAngiogenesisState2.eps
print -dpng tumorAntiAngiogenesisState2.png

figure(3);
pp = plot(solution.phase(1).time,solution.phase(1).control,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on;
print -depsc2 tumorAntiAngiogenesisControl.eps
print -dpng tumorAntiAngiogenesisControl.png

figure(4);
pp = plot(solution.phase(1).time,solution.phase(1).costate(:,1),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_p(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on;
print -depsc2 tumorAntiAngiogenesisCostate1.eps
print -dpng tumorAntiAngiogenesisConstate1.png

figure(5);
pp = plot(solution.phase(1).time,solution.phase(1).costate(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_q(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on;
print -depsc2 tumorAntiAngiogenesisCostate2.eps
print -dpng tumorAntiAngiogenesisConstate2.png

figure(6);
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
grid on;
print -depsc2 tumorAntiAngiogenesisMeshRefinement.eps
print -dpng tumorAntiAngiogenesisMeshRefinement.png

figure(7);
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
grid on;
print -depsc2 tumorAntiAngiogenesisMeshRefinementTime.eps
print -dpng tumorAntiAngiogenesisMeshRefinementTime.png

