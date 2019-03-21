function comparison_trig_paper()
close all
cd ./data
load('results_aoa_approx.mat')
data1 = out.setCONT(end).CONT(end).sol;
load('results_aoadot.mat')
data2 = out.setCONT(end).CONT(end).sol;
cd ..

re = 6378000;
H_scale = 7500;
rho0 = 1.2;
A = 557.4;
T = 2e6;
Isp = 400;
mu = 3.986e14;
g0 = 9.80665;
aoamax = 50;

% Approx Approach
tau_appr = data1.x;
time_appr = data1.parameters(1)*tau_appr;
alt_appr = data1.y(1,:)/1000;
h_appr = data1.y(1,:);
down_appr = data1.y(2,:)*6378;
v_appr = data1.y(3,:);
fpa_appr = data1.y(4,:)*180/pi;
gam_appr = data1.y(4,:);
mass_appr = data1.y(5,:);
aoa_appr = aoamax*sin(data1.control(1,:));
aoa_rad_appr = aoa_appr*pi/180;

lamalt_appr = data1.y(6,:);
lamdown_appr = data1.y(7,:);
lamv_appr = data1.y(8,:);
lamfpa_appr = data1.y(9,:);
lammass_appr = data1.y(10,:);

for count = 1:1:length(aoa_appr)
    if count == length(aoa_appr) || count == length(aoa_appr)-1
        aoadot_appr(1,count) = 0;
    else
        aoadot_appr(1,count) = (aoa_appr(1,count+1)-aoa_appr(1,count))/...
                               (time_appr(1,count+1)-time_appr(1,count));
    end
end

r_appr = re + h_appr;
Cl_appr = 0.4639*aoa_rad_appr-0.0278;
Cd_appr = 0.3216*aoa_rad_appr.^2-0.0305*aoa_rad_appr+0.03;
rho_appr = rho0*exp(-h_appr/H_scale);
q_appr = 0.5*rho_appr.*v_appr.^2;
L_appr = A.*q_appr.*Cl_appr;
D_appr = A.*q_appr.*Cd_appr;

ham_appr1 = lamalt_appr.*(v_appr.*sin(gam_appr));
ham_appr2 = lamdown_appr.*(v_appr.*cos(gam_appr)./r_appr);
ham_appr3 = lamv_appr.*((T-D_appr)./mass_appr-mu.*sin(gam_appr)./r_appr.^2);
ham_appr4 = lamfpa_appr.*((T.*(aoa_rad_appr)+L_appr)./(mass_appr.*v_appr)+(v_appr./r_appr - mu./(v_appr.*r_appr.^2)).*cos(gam_appr));
ham_appr5 = lammass_appr.*(-T./(g0.*Isp));
ham_appr = ham_appr1 + ham_appr2 + ham_appr3 + ham_appr4 + ham_appr5;

% Trig Approach
tau_trig = data2.x;
time_trig = data2.parameters(1)*tau_trig;
alt_trig = data2.y(1,:)/1000;
h_trig = data2.y(1,:);
down_trig = data2.y(2,:)*6378;
v_trig = data2.y(3,:);
fpa_trig = data2.y(4,:)*180/pi;
gam_trig = data2.y(4,:);
mass_trig = data2.y(5,:);
aoa_trig = data2.y(6,:)*180/pi;
aoa_rad_trig = data2.y(6,:);

lamalt_trig = data2.y(7,:);
lamdown_trig = data2.y(8,:);
lamv_trig = data2.y(9,:);
lamfpa_trig = data2.y(10,:);
lammass_trig = data2.y(11,:);
lamaoa_trig = data2.y(12,:);

aoadot_trig = 5*sin(data2.control(1,:));
aoadot_rad_trig = aoadot_trig.*pi/180;

r_trig = re + h_trig;
Cl_trig = 0.4639*aoa_rad_trig-0.0278;
Cd_trig = 0.3216*aoa_rad_trig.^2-0.0305*aoa_rad_trig+0.03;
rho_trig = rho0*exp(-h_trig/H_scale);
q_trig = 0.5*rho_trig.*v_trig.^2;
L_trig = A.*q_trig.*Cl_trig;
D_trig = A.*q_trig.*Cd_trig;

ham_trig1 = lamalt_trig.*(v_trig.*sin(gam_trig));
ham_trig2 = lamdown_trig.*(v_trig.*cos(gam_trig)./r_trig);
ham_trig3 = lamv_trig.*((T.*cos(aoa_rad_trig)-D_trig)./mass_trig-mu.*sin(gam_trig)./r_trig.^2+cos(data2.control(1,:)));
ham_trig4 = lamfpa_trig.*((T.*sin(aoa_rad_trig)+L_trig)./(mass_trig.*v_trig)+(v_trig./r_trig - mu./(v_trig.*r_trig.^2)).*cos(gam_trig));
ham_trig5 = lammass_trig.*(-T./(g0.*Isp));
ham_trig6 = lamaoa_trig.*aoadot_rad_trig;
ham_trig = ham_trig1 + ham_trig2 + ham_trig3 + ham_trig4 + ham_trig5 + ham_trig6;

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%
lw = 1.5;

% States Time History Plot
figure(1)
subplot(3,2,1)
h1 = plot(time_appr, alt_appr, 'b','markersize', 3, 'linewidth', lw);
ylabel('Altitude [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, alt_trig, 'r-.','markersize', 3, 'linewidth', lw);
% h3 = plot(time_appr, 50*ones(size(time_appr)), 'k--','markersize', 3, 'linewidth', 1);
xlim([0 146])
ylim([0 55])
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

subplot(3,2,2)
h1 = plot(time_appr, down_appr, 'b','markersize', 3, 'linewidth', lw);
ylabel('Downrange [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, down_trig, 'r-.','markersize', 3, 'linewidth', lw);
xlim([0 146])
ylim([0 350])
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

subplot(3,2,3)
h1 = plot(time_appr, v_appr/1000, 'b','markersize', 3, 'linewidth', lw);
ylabel('Velocity [km/s]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, v_trig/1000, 'r-.','markersize', 3, 'linewidth', lw);
xlim([0 146])
ylim([1 3.65])
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

subplot(3,2,4)
h1 = plot(time_appr, fpa_appr, 'b','markersize', 3, 'linewidth', lw);
ylabel('Flight Path Angle [deg]', 'fontSize', 16, 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16, 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, fpa_trig, 'r-.','markersize', 3, 'linewidth', lw);
xlim([0 146])
ylim([-74 22.75])
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

subplot(3,2,5)
h1 = plot(time_appr, mass_appr/1e5, 'b','markersize', 3, 'linewidth', lw);
ylabel('Mass X 10^{5} [kg]', 'fontSize', 16, 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16, 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, mass_trig/1e5, 'r-.','markersize', 3, 'linewidth', lw);
xlim([0 146])
ylim([0.6 1.4])
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

subplot(3,2,6)
h1 = plot(time_appr, aoa_appr, 'b','markersize', 3, 'linewidth', lw);
ylabel('Angle of Attack [deg]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, aoa_trig, 'r-.','markersize', 3, 'linewidth', lw);
xlim([0 146])
ylim([-50 2.4])
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

% figure(2)
% subplot(3,2,1)
% h1 = plot(time_appr, lamalt_appr*10000, 'b','markersize', 3, 'linewidth', lw);
% ylabel('\lambda_{h} X 10^{-4} [s/m]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_trig, lamalt_trig, 'r*','markersize', 3, 'linewidth', lw);
% set(gca,'FontSize',16,'FontWeight' , 'bold');
% legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)
% 
% subplot(3,2,2)
% h1 = plot(time_appr, lamdown_appr/1000, 'b','markersize', 3, 'linewidth', lw);
% ylabel('\lambda_{\theta} X 1000 [s/rad]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_trig, lamdown_trig/1000, 'r*','markersize', 3, 'linewidth', lw);
% set(gca,'FontSize',16,'FontWeight' , 'bold');
% legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)
% 
% subplot(3,2,3)
% h1 = plot(time_appr, lamv_appr, 'b','markersize', 3, 'linewidth', lw);
% ylabel('\lambda_{v} [s^2/m]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_trig, lamv_trig, 'r*','markersize', 3, 'linewidth', lw);
% set(gca,'FontSize',16,'FontWeight' , 'bold');
% legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)
% 
% subplot(3,2,4)
% h1 = plot(time_appr, lamfpa_appr, 'b','markersize', 3, 'linewidth', lw);
% ylabel('\lambda_{\gamma} [s/rad]', 'fontSize', 16, 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16, 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_trig, lamfpa_trig, 'r*','markersize', 3, 'linewidth', lw);
% set(gca,'FontSize',16,'FontWeight' , 'bold');
% legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)
% 
% subplot(3,2,5)
% h1 = plot(time_appr, lammass_appr*10000, 'b','markersize', 3, 'linewidth', lw);
% ylabel('\lambda_{m} X 10^{-4} [s/kg]', 'fontSize', 16, 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16, 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_trig, lammass_trig, 'r*','markersize', 3, 'linewidth', lw);
% set(gca,'FontSize',16,'FontWeight' , 'bold');
% legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)
% 
% subplot(3,2,6)
% h1 = plot(time_appr, aoadot_appr, 'b','markersize', 3, 'linewidth', lw);
% ylabel('\alpha Rate [deg/s]', 'fontSize', 16, 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16, 'fontWeight' , 'bold')
% hold on
% h2 = plot(time_trig, aoadot_trig, 'r*','markersize', 3, 'linewidth', lw);
% set(gca,'FontSize',16,'FontWeight' , 'bold');
% legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

figure(2)
subplot(1,2,1)
h1 = plot(time_appr, ham_appr, 'b','markersize', 3, 'linewidth', lw);
ylabel('Hamiltonian [nd]', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, ham_trig, 'r-.','markersize', 3, 'linewidth', lw);
xlim([0 146])
set(gca,'FontSize',20,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)

subplot(1,2,2)
h1 = plot(time_appr, (ham_appr+1)*1e8, 'b','markersize', 3, 'linewidth', lw);
ylabel('[Hamiltonian + 1] X 10^{-8} [nd]', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
hold on
h2 = plot(time_trig, (ham_trig+1)*1e8, 'r-.','markersize', 3, 'linewidth', lw);
xlim([0 146])
ylim([-250 460])
set(gca,'FontSize',20,'FontWeight' , 'bold');
legend([h1 h2],{'Approximation','Trigonometrization'},'fontSize', 14)
return