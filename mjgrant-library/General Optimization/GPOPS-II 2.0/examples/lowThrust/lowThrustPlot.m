%------------------------------%
% Extract Solution from Output %
%------------------------------%
time = solution.phase.time;
p = solution.phase.state(:,1);
f = solution.phase.state(:,2);
g = solution.phase.state(:,3);
h = solution.phase.state(:,4);
k = solution.phase.state(:,5);
L = solution.phase.state(:,6);
w = solution.phase.state(:,7);

[a,ecc,inc,Ome,ome,nu] = lowThrustMee2Coe(p,f,g,h,k,L);

ur = solution.phase.control(:,1);
ut = solution.phase.control(:,2);
uh = solution.phase.control(:,3);

tau = solution.parameter;

q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;

% ----------------------------------------------------------------------- %
% Trajectory in Cartesian
% ----------------------------------------------------------------------- % 
X = (r./s2).*(cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L));
Y = (r./s2).*(sin(L)-alpha2.*sin(L)+2.*h.*k.*cos(L));
Z = (2.*r./s2).*(h.*sin(L)-k.*cos(L));

%---------------%
% Plot Solution %
%---------------%
timehr = time/(60*60);

xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 .2 0]);

figure(1);
h1 = gca;
earth_sphere(h1,'ftmod'); hold on;
plot3(X./10000000,Y./10000000,Z./10000000,'r','LineWidth',1.25); 
axis([-6 6 -6 6 -3 10]); grid on;
xl = xlabel('$x(t)$~(ft$\times 10^{7}$)','Interpreter','LaTeX'); 
yl = ylabel('$y(t)$~(ft$\times 10^{7}$)','Interpreter','LaTeX'); 
zl = zlabel('$y(t)$~(ft$\times 10^{7}$)','Interpreter','LaTeX'); 
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(zl,'FontSize',16);
set(gca,'XTick',[-5 0 5]);
set(gca,'YTick',[-5 0 5]);
set(gca,'ZTick',[0 5 10]);
set(gca,'FontSize',16,'FontName','Times');
hold off;
print -dpng lowThrustOptimalTransfer.png;

figure(2);
plot(timehr,p./1000000,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(p/1000000) max(p/1000000)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$p(t)$~(ft$\times{10}^6$)','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[30 40 50]);
print -dpng lowThrustp.png;

figure(3);
plot(timehr,f,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(f) max(f)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$f(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[-0.1 -0.05 -0.0]);
print -dpng lowThrustf.png;

figure(4);
plot(timehr,g,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(g) max(g)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$g(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[0.0 0.2 0.4 0.6]);
print -dpng lowThrustg.png;

figure(5);
plot(timehr,h,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(h) max(h)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$h(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[-0.6 -0.5 -0.4 -0.3]);
print -dpng lowThrusth.png;

figure(6);
plot(timehr,k,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(k) max(k)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$k(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[-0.1 -0.05 0.0 0.05]);
print -dpng lowThrustk.png;

figure(7);
plot(timehr,L./(2*pi),'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(L/(2*pi)) max(L/(2*pi))]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$L(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[2 4 6 8]);
print -dpng lowThrustL.png;

figure(8);
plot(timehr,ur,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(ur) max(ur)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$u_r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[-0.5 0 0.5]);
print -dpng lowThrustur.png;

figure(9);
plot(timehr,ut,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(ut) max(ut)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$u_t(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[-0.5 0 0.5]);
print -dpng lowThrustut.png;

figure(10);
plot(timehr,uh,'-','LineWidth',1.25); grid on;
axis ([min(timehr) max(timehr) min(uh) max(uh)]);
xl = xlabel('$t$~(hr)','Interpreter','LaTeX');
yl = ylabel('$u_h(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'YTick',[-0.5 0 0.5]);
print -dpng lowThrustuh.png;
