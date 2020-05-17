% Comparison Plot File
clear;close;clc
cd data
load('results')
load('msledl_gpops_ac')

data_bvp4c = out.setCONT(end).CONT(end).sol;
tau_bvp4c = data_bvp4c.x;
tf_bvp4c = data_bvp4c.parameters(1);
time_bvp4c = tau_bvp4c.*tf_bvp4c;
h_bvp4c = data_bvp4c.y(1,:);
v_bvp4c = data_bvp4c.y(2,:);
gam_bvp4c = data_bvp4c.y(3,:);
londot = v_bvp4c.*cos(gam_bvp4c);
lon = trapz(time_bvp4c,londot);
downrange_bvp4c = lon.*time_bvp4c./max(time_bvp4c);

lamh_bvp4c = data_bvp4c.y(4,:);
lamv_bvp4c = data_bvp4c.y(5,:);
lamgam_bvp4c = data_bvp4c.y(6,:);
control_bvp4c = data_bvp4c.control;
bank_bvp4c = acosd(0.683*sin(control_bvp4c) + 0.183);

H = 9354;
rho0 = 0.0158;
rn = 0.6;
k = 1.9027e-4;
Cd = 1.45;
LbyD = 0.24;
Cl = Cd.*LbyD;
mass = 3300;
g = 9.81;
Aref = 15.9;
epsQ = 1e-6;
epsc = 1e-6;
epsq = 1e-6;
epsg = 1e-6;
rm = 3397000;
muu = 4.284*1e13;
Qmax = 70e4;
gmax = 5;
qmax = 10e3;
umin = -0.5;
umax = 0.866;

q_bvp4c = 0.5*(rho0.*exp(-h_bvp4c./H)).*v_bvp4c.^2;
Q_bvp4c = k.*sqrt((rho0.*exp(-h_bvp4c./H))./rn).*v_bvp4c.^3;
g_bvp4c = sqrt((1/2.*(rho0.*exp(-h_bvp4c./H)).*v_bvp4c.^2*Cd.*LbyD.*Aref).^2 + (1/2.*(rho0.*exp(-h_bvp4c./H)).*v_bvp4c.^2.*Cd.*Aref).^2)./(mass.*g);
h_q_bvp4c = -H.*log((2*9000)./(rho0*v_bvp4c.^2));
h_Q_bvp4c = -H.*log((rn.*700000^2)./(rho0.*k.^2.*v_bvp4c.^6));
h_g_bvp4c = -H.*log((10.*mass.*9.81)./(v_bvp4c.^2.*Aref.*sqrt(Cl.^2+Cd.^2)));

ham_bvp4c = epsQ./cos((0.5*k*v_bvp4c.^3.*pi.*((rho0.*exp(-h_bvp4c./H))./rn).^(1/2))./Qmax)...
            - lamv_bvp4c.*((muu.*sin(gam_bvp4c))./(h_bvp4c + rm).^2 ...
            + (Aref*Cd*rho0*v_bvp4c.^2.*exp(-h_bvp4c./H))/(2*mass))...
            + epsq./cos((0.25*rho0*v_bvp4c.^2.*pi.*exp(-h_bvp4c./H))./qmax)...
            + epsg./cos((0.5*pi*((Aref^2*Cd^2*rho0^2*v_bvp4c.^4.*exp(-(2.*h_bvp4c)./H))./4 ...
            + (Aref^2*Cd^2*LbyD^2*rho0^2*v_bvp4c.^4.*exp(-(2*h_bvp4c)./H))./4).^(1/2))./(g*gmax*mass))...
            + lamgam_bvp4c.*((v_bvp4c.*cos(gam_bvp4c))./(h_bvp4c + rm)...
            - (muu*cos(gam_bvp4c))./(v_bvp4c.*(h_bvp4c + rm).^2)...
            + (Aref*Cd*LbyD*rho0*v_bvp4c.*exp(-h_bvp4c./H).*(0.683.*sin(control_bvp4c) + 0.183))./(2*mass))...
            + epsc*cos(control_bvp4c) + lamh_bvp4c.*v_bvp4c.*sin(gam_bvp4c);

sf_bvp4c = 0.5*Aref*Cl*rho0*v_bvp4c.*exp(-h_bvp4c./H).*lamgam_bvp4c./mass;

% GPOPS-II
time_gpops = time;
h_gpops = h;
downrange_gpops = Downrange;
v_gpops = v;
gam_gpops = gam*pi/180;
bank_gpops = bank;
q_gpops = Q;
Q_gpops = q;
g_gpops = gload;
h_q_gpops = h_q;
h_Q_gpops = h_heat;
h_g_gpops = h_g;

lamh_gpops = solution.phase.costate(:,1);
lamv_gpops = solution.phase.costate(:,2);
lamgam_gpops = solution.phase.costate(:,3);
cd ..

qmax = ones(size(time_bvp4c))*10;
Qmax = ones(size(time_bvp4c))*70;
gmax = ones(size(time_bvp4c))*5;

ham_gpops = - lamv_gpops.*((muu.*sin(gam_gpops))./(h_gpops + rm).^2 ...
            + (Aref*Cd*rho0*v_gpops.^2.*exp(-h_gpops./H))/(2*mass))...
            + lamgam_gpops.*((v_gpops.*cos(gam_gpops))./(h_gpops + rm)...
            - (muu*cos(gam_gpops))./(v_gpops.*(h_gpops + rm).^2)...
            + (Aref*Cd*LbyD*rho0*v_gpops.*exp(-h_gpops./H).*cosd(bank_gpops))./(2*mass))...
            + lamh_gpops.*v_gpops.*sin(gam_gpops);

sf_gpops = 0.5*Aref*Cl*rho0*v_gpops.*exp(-h_gpops./H).*lamgam_gpops./mass;

% %---------------%
% % Plot Solution %
% %---------------%
figure(1)
h1 = plot(v_bvp4c/1000,h_bvp4c/1000,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Velocity (km/s)');
yl = ylabel('Altitude (km)');
set(xl,'FontSize',20,'FontWeight' , 'bold');
set(yl,'FontSize',20,'FontWeight' , 'bold');
xlim([0.5 6.5])
ylim([0 130])
hold on 
h2 = plot(v_gpops/1000,h_gpops/1000,'rx', 'markersize', 3, 'linewidth', 2);
h3 = plot(v_gpops/1000,h_q_gpops,'color',[0 0.5 0],'LineStyle',':','markersize', 1, 'linewidth', 2);
h4 = plot(v_gpops/1000,h_Q_gpops,'k-+', 'markersize', 1, 'linewidth', 2);
h5 = plot(v_gpops/1000,h_g_gpops,'color',[.61 .51 .74],'LineStyle','--', 'markersize', 1, 'linewidth', 2);
legend([h1 h2 h3 h4 h5],{'UTM','PSM','Dynamic Pressure Constraint','Heat-Rate Constraint','g-Load Constraint'},'fontSize', 12)
set(gca,'FontSize',20,'FontWeight' , 'bold');

figure(2)
h1 = plot(downrange_bvp4c/1000,h_bvp4c/1000,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Downrange (km)');
yl = ylabel('Altitude (km)');
set(xl,'FontSize',20,'FontWeight' , 'bold');
set(yl,'FontSize',20,'FontWeight' , 'bold');
xlim([0 1100])
ylim([0 130])
hold on 
h2 = plot(downrange_gpops,h_gpops/1000,'rx', 'markersize', 3, 'linewidth', 2);
h3 = plot(downrange_gpops,h_q_gpops,'color',[0 0.5 0],'LineStyle',':', 'markersize', 1, 'linewidth', 2);
h4 = plot(downrange_gpops,h_Q_gpops,'k-+', 'markersize', 1, 'linewidth', 2);
h5 = plot(downrange_gpops,h_g_gpops,'color',[.61 .51 .74],'LineStyle','--', 'markersize', 1, 'linewidth', 2);
legend([h1 h2 h3 h4 h5],{'UTM','PSM','Dynamic Pressure Constraint','Heat-Rate Constraint','g-Load Constraint'},'fontSize', 12)
set(gca,'FontSize',20,'FontWeight' , 'bold');

figure(3)
subplot(2,2,1)
h1 = plot(time_bvp4c,bank_bvp4c,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('Bank Angle (deg)');
set(xl,'FontSize',20,'FontWeight' , 'bold');
set(yl,'FontSize',20,'FontWeight' , 'bold');
xlim([0 320])
ylim([25 125])
hold on 
h2 = plot(time_gpops,bank_gpops,'r-x', 'markersize', 3, 'linewidth', 1);
legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

subplot(2,2,2)
h1 = plot(time_bvp4c,q_bvp4c/1e3,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('Dynamic Pressure (kPa)');
set(xl,'FontSize',20,'FontWeight' , 'bold');
set(yl,'FontSize',20,'FontWeight' , 'bold');
xlim([0 320])
ylim([0 10.1])
hold on 
h2 = plot(time_gpops,q_gpops/1e3,'rx', 'markersize', 3, 'linewidth', 2);
h3 = plot(time_bvp4c,qmax,'k--', 'markersize', 1, 'linewidth', 2);
legend([h1 h2 h3],{'UTM','PSM','q_{MAX}'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

subplot(2,2,3)
h1 = plot(time_bvp4c,Q_bvp4c/1e4,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('Heat Rate (W/cm^{2})');
set(xl,'FontSize',20,'FontWeight' , 'bold');
set(yl,'FontSize',20,'FontWeight' , 'bold');
xlim([0 320])
ylim([0 70.1])
hold on 
h2 = plot(time_gpops,Q_gpops,'rx', 'markersize', 3, 'linewidth', 2);
h3 = plot(time_bvp4c,Qmax,'k--', 'markersize', 1, 'linewidth', 2);
legend([h1 h2 h3],{'UTM','PSM','Q-Rate_{MAX}'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

subplot(2,2,4)
h1 = plot(time_bvp4c,g_bvp4c,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('G-Load');
set(xl,'FontSize',20,'FontWeight' , 'bold');
set(yl,'FontSize',20,'FontWeight' , 'bold');
xlim([0 320])
ylim([0 5.1])
hold on 
h2 = plot(time_gpops,g_gpops,'rx', 'markersize', 3, 'linewidth', 2);
h3 = plot(time_bvp4c,gmax,'k--', 'markersize', 1, 'linewidth', 2);
legend([h1 h2 h3],{'UTM','PSM','g_{MAX}'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

figure(4)
subplot(1,2,1)
h1 = plot(time_bvp4c,ham_bvp4c,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('Hamiltonian (m/s)');
set(xl,'FontSize',18,'FontWeight' , 'bold');
set(yl,'FontSize',18,'FontWeight' , 'bold');
xlim([0 50])
hold on 
h2 = plot(time_gpops,ham_gpops,'r-x', 'markersize', 3, 'linewidth', 2);
legend([h1 h2],{'UTM','PSM'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(1,2,2)
h1 = plot(time_bvp4c,sf_bvp4c*1e3,'b', 'markersize', 2, 'linewidth', 2);
xl = xlabel('Time (s)');
yl = ylabel('Switching Function (m/(s deg)) X 10^{-3}');
set(xl,'FontSize',18,'FontWeight' , 'bold');
set(yl,'FontSize',18,'FontWeight' , 'bold');
xlim([0 50])
hold on 
h2 = plot(time_gpops,sf_gpops*1e3,'r-x', 'markersize', 3, 'linewidth', 2);
legend([h1 h2],{'UTM','PSM'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

% figure(5)
% h1 = plot(time_bvp4c,v_bvp4c/1000,'b', 'markersize', 2, 'linewidth', 2);
% xl = xlabel('Time (s)');
% yl = ylabel('Velocity (km/s)');
% set(xl,'FontSize',20,'FontWeight' , 'bold');
% set(yl,'FontSize',20,'FontWeight' , 'bold');
% hold on 
% h2 = plot(time_gpops,v_gpops/1000,'r-x', 'markersize', 3, 'linewidth', 2);
% legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
% set(gca,'FontSize',20,'FontWeight' , 'bold');

