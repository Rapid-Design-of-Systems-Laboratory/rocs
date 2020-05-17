function bounded_control_comp_plots()
clc;close all; clear
H = 11100;
rho0 = 0.02;
Aref = 250;
mass = 92080;
mu = 42828.371901.*1e9;
re = 3397000;
alfamax = 20*pi/180;

Cl1     = 1.6756;
Cl0     = -0.2070;
Cd2     = 2.04;
Cd1     = -0.3529;
Cd0     = 0.0785;

cd ./data
% load('results_bvp4c.mat')
load('results.mat')
data1 = out.setCONT(end).CONT(end).sol;
% UTM Result
h_bvp4c = data1.y(1,:)./1000;
down_bvp4c = (data1.y(2,:).*3397);
v_bvp4c = data1.y(3,:);
gam_bvp4c = data1.y(4,:).*180./pi;
u_bvp4c = 20*sin(data1.control(1,:));
tau_bvp4c = data1.x;
% data1.parameters(1) = 1;
time_bvp4c = (data1.parameters(1).*tau_bvp4c); % Scale unbounded time by 
lamh_bvp4c = data1.y(5,:);
lamthetta_bvp4c = data1.y(6,:);
lamv_bvp4c = data1.y(7,:);
lamgam_bvp4c = data1.y(8,:);
h = data1.y(1,:);
v = data1.y(3,:);
gam = data1.y(4,:);
lamH = data1.y(5,:);
lamTHETTA = data1.y(6,:);
lamV = data1.y(7,:);
lamGAM = data1.y(8,:);
alfatrig = data1.control(1,:);
ham_bvp4c = (lamGAM.*(cos(gam).*(v./(h + re) - mu./(v.*(h + re).^2)) + (Aref.*rho0.*v.*exp(-h./H).*(Cl0 + Cl1.*alfamax.*sin(alfatrig)))./(2.*mass)) - lamV.*((mu.*sin(gam))./(h + re).^2 + (Aref.*rho0.*v.^2.*exp(-h./H).*(Cd0 + Cd2.*alfamax.^2.*sin(alfatrig).^2 + Cd1.*alfamax.*sin(alfatrig)))./(2.*mass)) + lamH.*v.*sin(gam) + (lamTHETTA.*v.*cos(gam))./(h + re));

load('results_gpops.mat')
% PSM Result
h_gpops = h;
down_gpops = thetta*3397;
v_gpops = v;
gam_gpops = gam;
u_gpops = aoa;
time_gpops = time;
lamh_gpops = lamH;
lamthetta_gpops = lamTHETTA;
lamv_gpops = lamV;
lamgam_gpops = lamGAM; 
cd ..

alfalb = -20*ones(size(time_gpops));
alfaub = 20*ones(size(time_gpops));
%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

figure(1)
subplot(2,2,1)
h1 = plot(v_bvp4c./1000,h_bvp4c, 'b','markersize', 2, 'linewidth', 3);
ylabel('Altitude [km]', 'fontSize', 20 , 'fontWeight' , 'bold')
xlabel('Velocity [km/s]', 'fontSize', 20 , 'fontWeight' , 'bold')
hold on
h2 = plot(v_gpops./1000,h_gpops, 'rx','markersize', 7, 'linewidth', 1);
legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

subplot(2,2,2)
h1 = plot(down_bvp4c,h_bvp4c, 'b','markersize', 2, 'linewidth', 3);
ylabel('Altitude [km]', 'fontSize', 20 , 'fontWeight' , 'bold')
xlabel('Downrange [km]', 'fontSize', 20 , 'fontWeight' , 'bold')
hold on
h2 = plot(down_gpops,h_gpops, 'rx','markersize', 7, 'linewidth', 1);
legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

subplot(2,2,3)
h1 = plot(time_bvp4c,gam_bvp4c, 'b','markersize', 2, 'linewidth', 3);
ylabel('Flight Path Angle [deg]', 'fontSize', 20 , 'fontWeight' , 'bold')
xlabel('Time', 'fontSize', 20 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_gpops,gam_gpops, 'rx','markersize', 7, 'linewidth', 1);
legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

subplot(2,2,4)
h1 = plot(time_bvp4c,u_bvp4c, 'b','markersize', 2, 'linewidth', 3);
ylabel('Angle of Attack [deg]', 'fontSize', 20 , 'fontWeight' , 'bold')
xlabel('Time', 'fontSize', 20 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_gpops,u_gpops, 'rx','markersize', 7, 'linewidth', 1);
h4 = plot(time_gpops,alfaub, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h1 h2 h4],{'UTM','PSM','Upper Bound'},'fontSize', 16)
set(gca,'FontSize',20,'FontWeight' , 'bold');

% figure(2)
% subplot(2,2,1)
% h1 = plot(time_bvp4c,lamh_bvp4c.*1000, 'b','markersize', 3, 'linewidth', 2);
% ylabel('\lambda_{h} [s/km]', 'fontSize', 20 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_gpops,lamh_gpops.*1000, 'r-+','markersize', 3, 'linewidth', 1);
% legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
% set(gca,'FontSize',20,'FontWeight' , 'bold');
% 
% subplot(2,2,2)
% h1 = plot(time_bvp4c,lamthetta_bvp4c*1e6, 'b','markersize', 1, 'linewidth', 2);
% ylabel('\lambda_{\theta} [s/rad] X 10^{-6}', 'fontSize', 20 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_gpops,lamthetta_gpops*1e6, 'r-+','markersize', 3, 'linewidth', 1);
% legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
% set(gca,'FontSize',20,'FontWeight' , 'bold');
% 
% subplot(2,2,3)
% h1 = plot(time_bvp4c,lamv_bvp4c, 'b','markersize', 3, 'linewidth', 1.5);
% ylabel('\lambda_{v} [s^2/m]', 'fontSize', 20 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_gpops,lamv_gpops, 'r-+','markersize', 3, 'linewidth', 1.5);
% legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
% set(gca,'FontSize',20,'FontWeight' , 'bold');
% 
% subplot(2,2,4)
% h1 = plot(time_bvp4c,lamgam_bvp4c./1000, 'b','markersize', 1, 'linewidth', 1.5);
% ylabel('\lambda_{\gamma} [s/rad] X 10^{3}', 'fontSize', 20 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_gpops,lamgam_gpops./1000, 'r-+','markersize', 1, 'linewidth', 1.5);
% legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
% set(gca,'FontSize',20,'FontWeight' , 'bold');
% 
% figure(3)
% h1 = plot(time_bvp4c,ham_bvp4c, 'b','markersize', 3, 'linewidth', 1.5);
% ylabel('Hamiltonian [nd]', 'fontSize', 20 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_gpops,ham_gpops, 'r-+','markersize', 3, 'linewidth', 1.5);
% legend([h1 h2],{'UTM','PSM'},'fontSize', 16)
% set(gca,'FontSize',20,'FontWeight' , 'bold');
return