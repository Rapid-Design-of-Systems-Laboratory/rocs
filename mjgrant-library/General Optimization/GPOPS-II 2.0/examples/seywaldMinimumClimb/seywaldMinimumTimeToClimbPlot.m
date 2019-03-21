%-------------------------------------------------------------------------%
%                             Extract Solution                            %
%-------------------------------------------------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
state = solution.phase(1).state;
control = solution.phase(1).control;
for i=1:length(output.meshhistory);
  mesh(i).points = [0 cumsum(output.meshhistory(i).result.setup.mesh.phase.fraction)];
  mesh(i).iteration = i*ones(size(mesh(i).points));
end;

%-------------------------------------------------------------------------%
%                              Plot Solution                              %
%-------------------------------------------------------------------------%
figure
pp = plot(time,state(:,1)/1000,'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$h(t)$~(km)','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
% print -depsc2 moonLanderAltitude.eps
% print -dpng moonLanderAltitude.png

figure
pp = plot(time,state(:,2),'-o');
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$v(t)$ (m$\cdot$s${}^{-1})$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
% print -depsc2 moonLanderVelocity.eps
% print -dpng moonLanderVelocity.png

figure
pp = plot(state(:,2),state(:,1)/1000,'-o');
xl = xlabel('$v$ (m$\cdot$s${}^{-1})$','Interpreter','LaTeX');
yl = ylabel('$h$~(km)','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
% print -depsc2 moonLanderVelocity.eps
% print -dpng moonLanderVelocity.png

figure
pp = plot(time,control,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
% print -depsc2 moonLanderControl.eps
% print -dpng moonLanderControl.png

figure
tf = time(end);
for i=1:length(mesh);
  pp = plot(mesh(i).points*tf,mesh(i).iteration,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Mesh Point Location (Fraction of Interval)');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16);
grid on;
% print -depsc2 moonLanderMeshRefinement.eps
% print -dpng moonLanderMeshRefinement.png
