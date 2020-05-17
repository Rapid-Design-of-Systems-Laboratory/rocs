function comparison_plots()
close all; clc
cd ./data

load('results.mat')
data_u = out.setCONT(end).CONT(end).sol;
for i = 1:1:length(out.setCONT(end).CONT)
    data = out.setCONT(end).CONT(i);
    obj_u(i) = data.sol.y(3,end);
    error_u(i) = data.in.const.epsilonh{1,1};
end
cd ..

rho0     = 1.2; % Sea-level atmospheric density
H        = 7500; % Scale height for atmosphere of Earth
qmax     = 2000e4; 
k        = 1.74153e-4; % heat rate coefficient
rn       = 1/12*0.3048; % Nose radius

% UTM
h_u = data_u.y(1,:)/1000;
down_u = 6378*data_u.y(2,:);
v_u = data_u.y(3,:);
% gam_u = data_u.y(4,:)*180/pi;
% lamh_u = data_u.y(5,:);
% lamthetta_u = data_u.y(6,:);
% lamv_u = data_u.y(7,:);
% lamgam_u = data_u.y(8,:);
alfa_u = 40*sin(data_u.control(1,:));
tau_u = data_u.x;
tf_u = data_u.parameters(1);
time_u = tf_u*tau_u;
heat_u = k.*sqrt((rho0*exp(-h_u*1000/H))/rn).*v_u.^3*1e-4;
h_heat_u = -H*log(rn*qmax^2./(k^2*rho0*v_u.^6))/1000;
u_bound = 40*ones(size(time_u));
heat_bound = 2000*ones(size(time_u));
Cl_u = 1.5658*alfa_u*pi/180;
Cd_u = 1.6537*(alfa_u*pi/180).^2 + 0.0612;
LbyD_u = Cl_u./Cd_u;

% GPOPS-II
cd ./data
load('mc_gpops2.mat')

h_g = h;
down_g = downrange;
v_g = v;
gam_g = gam;
alfa_g = aoa;
time_g = time;
heat_g = k.*sqrt((rho0*exp(-h_g*1000/H))/rn).*v_g.^3*1e-4;
% h_heat_g = -H*log(rn*qmax^2./(k^2*rho0*v_g.^6))/1000;
Cl_g = 1.5658*alfa_g*pi/180;
Cd_g = 1.6537*(alfa_g*pi/180).^2 + 0.0612;
LbyD_g = Cl_g./Cd_g;
cd ..

% load('mc_results.mat')
error_u = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5];
obj_u = [2175.167506207363,2192.117215669861,2192.485730830485,2199.220695009121,2199.696751286896,2199.854022922838];

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

% States Time History Plot
figure(1)
subplot(1,2,1)
h1 = plot(v_u/1000,h_u, 'b*','markersize', 5, 'linewidth', 2);
ylabel('Altitude [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Velocity [km/s]', 'fontSize', 16 , 'fontWeight' , 'bold')
ylim([0 80])
xlim([2 5.1])
hold on
set(gca,'FontSize',18,'FontWeight' , 'bold');
h2 = plot(v_g/1000,h_g, 'r','markersize', 3, 'linewidth', 2);
h3 = plot(v_u/1000,h_heat_u, 'k--','markersize', 3, 'linewidth', 2);
legend([h1 h2 h3],{'UTM','PSM','Heat Rate Constraint'},'fontSize',14)

subplot(1,2,2)
semilogx(error_u,obj_u, 'b-*','markersize', 3, 'linewidth', 2, 'MarkerEdgeColor','red',...
    'MarkerFaceColor','red');
% plot(error_u,obj_u, 'b-*','markersize', 3, 'linewidth', 2);
ylabel('Objective Funtional, J [m/s]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('\epsilon [m/s^{2}]', 'fontSize', 16 , 'fontWeight' , 'bold')
%xlim([1e-5 1])
set(gca,'FontSize',18,'FontWeight' , 'bold');

% subplot(2,2,2)
% h1 = plot(time_u,v_u/1000, 'b-*','markersize', 3, 'linewidth', 2);
% ylabel('Velocity [km/s]', 'fontSize', 12 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
% hold on
% set(gca,'FontSize',12,'FontWeight' , 'bold');
% h2 = plot(time_g,v_g/1000, 'r--','markersize', 3, 'linewidth', 2);
% legend([h1 h2],{'UTM','GPOPS'},'fontSize',14)

% subplot(2,2,3)
% h1 = plot(time_u,gam_u, 'b-*','markersize', 3, 'linewidth', 2);
% ylabel('\gamma [deg]', 'fontSize', 12 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
% hold on
% set(gca,'FontSize',12,'FontWeight' , 'bold');
% h2 = plot(time_g,gam_g, 'r--','markersize', 3, 'linewidth', 2);
% legend([h1 h2],{'UTM','GPOPS'},'fontSize',14)

figure(2)
subplot(1,2,1)
h1 = plot(time_u,alfa_u, 'b*','markersize', 5, 'linewidth', 2);
ylabel('\alpha [deg]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
ylim([-2.5 41])
xlim([0 43])
hold on
set(gca,'FontSize',18,'FontWeight' , 'bold');
h2 = plot(time_g,alfa_g, 'r','markersize', 3, 'linewidth', 2);
h3 = plot(time_u,u_bound, 'k:','markersize', 3, 'linewidth', 2);
legend([h1 h2 h3],{'UTM','PSM','Upper Bound'},'fontSize',14)

subplot(1,2,2)
h1 = plot(time_u,heat_u, 'b*','markersize', 5, 'linewidth', 2);
ylabel('Heat Rate [W/cm^{2}]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
ylim([0 2100])
xlim([0 43])
hold on
set(gca,'FontSize',18,'FontWeight' , 'bold');
h2 = plot(time_g,heat_g, 'r','markersize', 3, 'linewidth', 2);
h3 = plot(time_u,heat_bound, 'k--','markersize', 3, 'linewidth', 2);
legend([h1 h2 h3],{'UTM','PSM','Upper Bound'},'fontSize',14)

figure(3)
h1 = plot(time_u,LbyD_u, 'b*','markersize', 5, 'linewidth', 2);
ylabel('L/D', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% ylim([-2.5 41])
xlim([0 43])
hold on
set(gca,'FontSize',18,'FontWeight' , 'bold');
h2 = plot(time_g,LbyD_g, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'UTM','PSM'},'fontSize',14)

% % Co-states Time History Plot
% figure(2)
% 
% subplot(2,2,1)
% plot(time_b,lamh_b*1000, 'b','markersize', 3, 'linewidth', 2)
% title('\lambda_{h} History', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('\lambda_{h} [s/km]', 'fontSize', 12 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
% grid on
% set(gca,'FontSize',12,'FontWeight' , 'bold');
% 
% subplot(2,2,2)
% plot(time_b,lamthetta_b, 'b','markersize', 3, 'linewidth', 2)
% title('\lambda_{\theta} History', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('\lambda_{\theta} [s/rad]', 'fontSize', 12 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
% grid on
% set(gca,'FontSize',12,'FontWeight' , 'bold');
% 
% subplot(2,2,3)
% plot(time_b,lamv_b, 'b','markersize', 3, 'linewidth', 2)
% title('\lambda_{\theta} History', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('\lambda_{v} [s^2/m]', 'fontSize', 12 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
% grid on
% set(gca,'FontSize',12,'FontWeight' , 'bold');
% 
% subplot(2,2,4)
% plot(time_b,lamgam_b, 'b','markersize', 3, 'linewidth', 2)
% title('\lambda_{\gamma} History', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('\lambda_{\gamma} [s/rad]', 'fontSize', 12 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
% grid on
% set(gca,'FontSize',12,'FontWeight' , 'bold');
