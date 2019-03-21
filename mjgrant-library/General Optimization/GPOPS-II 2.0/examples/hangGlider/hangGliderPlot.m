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

x  = solution.phase.state(:,1);
y  = solution.phase.state(:,2);
vx = solution.phase.state(:,3);
vy = solution.phase.state(:,4);
u  = solution.phase.control(:,1);

figure(1)
pp = plot(time,x,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$x(t)$~(m)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 HangGliderX.eps
print -dpng HangGliderX.png

figure(2)
pp = plot(time,y,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$y(t)$~(m)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 HangGliderY.eps
print -dpng HangGliderY.png

figure(3)
pp = plot(time,vx,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$v_x(t)$~(m$\cdot$s${}^{-1}$)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
grid on
print -depsc2 HangGliderVX.eps
print -dpng HangGliderVX.png

figure(4)
pp = plot(time,vy,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$v_y(t)$~(m$\cdot$s${}^{-1}$)','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 HangGliderVY.eps
print -dpng HangGliderVY.png

figure(5)
pp = plot(time,u,'-o', 'markersize', 7, 'linewidth', 1.5);
xl = xlabel('$t$~(s)','Interpreter','LaTeX');
yl = ylabel('$C_L(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 HangGliderU.eps
print -dpng HangGliderU.png

figure(6);
for i=1:length(mesh);
  pp = plot(mesh(i).meshPoints,mesh(i).iteration,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Mesh Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'FontSize',16);
axis([0,1,0,length(mesh)+1]);
grid on;

figure(7);
for i=1:length(mesh);
  pp = plot(mesh(i).time/mesh(i).time(end),mesh(i).iterationTime,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Collocation Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'FontSize',16);
axis([0,1,0,length(mesh)+1]);
grid on;
