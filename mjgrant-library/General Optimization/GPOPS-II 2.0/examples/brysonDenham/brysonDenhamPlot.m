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
%                         Compute Exact Solution                          %
%-------------------------------------------------------------------------%
index1 = find(time<=0.30000000001);
temp1 = find(time<=0.7);
index2 = (index1(end)+1:temp1(end)).';
index3 = (index2(end)+1:length(time)).';
t{1} = time(index1);
t{2} = time(index2);
t{3} = time(index3);

xExact{1}(:,1) = L*(1-(1-t{1}/(3*L)).^3);
xExact{2}(:,1) = L*ones(size(t{2}));
xExact{3}(:,1) = L*(1-(1-(1-t{3})/(3*L)).^3);
xExact{1}(:,2) = (1-t{1}/(3*L)).^2;
xExact{2}(:,2) = zeros(size(t{2}));
xExact{3}(:,2) = -(1-(1-t{3})/(3*L)).^2;
uExact{1}(:,1) = -(2/(3*L))*(1-t{1}/(3*L));
uExact{2}(:,1) = zeros(size(t{2}));
uExact{3}(:,1) = -(2/(3*L))*(1-(1-t{3})/(3*L));

tstar = [t{1}; t{2}; t{3}];
xstar = [xExact{1}; xExact{2}; xExact{3}];
ustar = [uExact{1}; uExact{2}; uExact{3}];

%-------------------------------------------------------------------------%
%                              Plot Solution                              %
%-------------------------------------------------------------------------%
figure(1);
pp = plot(time,state,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(x(t),v(t))$','Interpreter','LaTeX');
ll = legend('$x(t)$','$v(t)$');
set(ll,'Fontsize',18,'Interpreter','LaTeX');
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng brysonDenhamState.png

figure(2);
pp = plot(time,state(:,1),'-o',tstar,xstar(:,1),'-d');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$x(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng brysonDenhamX.png

figure(3);
pp = plot(time,state(:,2),'-o',tstar,xstar(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng brysonDenhamV.png

figure(4);
% pp = plot(time,control,'-o',tstar,ustar,'-d');
pp = plot(time,control,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
ll = legend('$u(t)$');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'Fontsize',18,'Interpreter','LaTeX');
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng brysonDenhamControl.png

figure(5);
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
set(gca,'YTick',0:length(mesh),'FontSize',16,'FontName','Times');
axis([0 tf 1 mesh(end).iteration(1)]);
grid on;
print -dpng brysonDenhamMeshHistory.png

figure(6);
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
grid on;
print -dpng brysonDenhamCollocHistory.png

