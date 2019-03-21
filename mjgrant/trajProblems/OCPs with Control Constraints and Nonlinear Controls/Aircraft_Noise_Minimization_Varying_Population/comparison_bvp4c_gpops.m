function comparison_bvp4c_gpops()
close all
%Constants
g  = 9.80665;
mass  = 7180/9.80665; % 7180
C1  = 0.226;
C2  = 5.2e6;
L = mass*g;

cd ./data
load('results_pop_xy.mat')
data1 = out.setCONT(end).CONT(end).sol;

% Traditional Approach
down_bvp4c = data1.y(1,:)/1000;
cross_bvp4c = data1.y(2,:)/1000;
alt_bvp4c = data1.y(3,:)/1000;
z_bvp4c = data1.y(3,:);
v_bvp4c = data1.y(4,:);
psi_bvp4c = data1.y(5,:)*180/pi;
psi_rad_bvp4c = data1.y(5,:);
gam_bvp4c = data1.y(6,:)*180/pi;
gam_rad_bvp4c = data1.y(6,:);
tau_bvp4c = data1.x;
time_bvp4c = data1.parameters(1)*tau_bvp4c;
lamX_bvp4c = data1.y(7,:);
lamY_bvp4c = data1.y(8,:);
lamZ_bvp4c = data1.y(9,:);
lamV_bvp4c = data1.y(10,:);
lamPSII_bvp4c = data1.y(11,:);
lamGAM_bvp4c = data1.y(12,:);
bank_bvp4c = 60*sin(data1.control(1,:));
bank_rad_bvp4c = bank_bvp4c.*pi/180;
aoa_bvp4c = 15*sin(data1.control(2,:));
aoa_rad_bvp4c = aoa_bvp4c.*pi/180;
T_bvp4c = 1560*sin(data1.control(3,:))+1860;
pop = 4600./(1+cross_bvp4c*1000);
Noise_bvp4c = 10*log10(((pop.*18.73.*T_bvp4c.^(5.2).*cos(gam_bvp4c.*pi./180))./(v_bvp4c.*(alt_bvp4c.*1000+50).^(2.5))));

D_bvp4c = C1.*v_bvp4c.^2 + C2./(v_bvp4c.^2);
Ft_bvp4c = T_bvp4c.*cos(aoa_rad_bvp4c) - D_bvp4c; % force along velocity vector
Fn_bvp4c = T_bvp4c.*sin(aoa_rad_bvp4c) + L; % force perpendicular to velocity vector
 
% Hamiltonian calculation
ham_bvp4c1 = lamX_bvp4c.*v_bvp4c.*cos(gam_rad_bvp4c).*cos(psi_rad_bvp4c); 
ham_bvp4c2 = lamY_bvp4c.*v_bvp4c.*cos(gam_rad_bvp4c).*sin(psi_rad_bvp4c);
ham_bvp4c3 = lamZ_bvp4c.*v_bvp4c.*sin(gam_rad_bvp4c);
ham_bvp4c4 = lamV_bvp4c.*(Ft_bvp4c./mass - g.*sin(gam_rad_bvp4c));
ham_bvp4c5 = lamPSII_bvp4c.*(Fn_bvp4c.*sin(bank_rad_bvp4c)./(mass.*cos(gam_rad_bvp4c).*v_bvp4c));
ham_bvp4c6 = lamGAM_bvp4c.*(Fn_bvp4c.*cos(bank_rad_bvp4c)./(mass.*v_bvp4c) - g.*cos(gam_rad_bvp4c)./v_bvp4c);
ham_bvp4c7 = T_bvp4c.^(5.2).*cos(gam_rad_bvp4c)./(v_bvp4c.*(z_bvp4c+50).^(2.5));
ham_bvp4c = ham_bvp4c1 + ham_bvp4c2 + ham_bvp4c3 + ham_bvp4c4 + ham_bvp4c5 + ham_bvp4c6 + ham_bvp4c7;

load('results_gpops_popy.mat')
down_gpops = down;
cross_gpops = crossrange;
alt_gpops = alt;
z_gpops = alt.*1000;
v_gpops = v;
psi_gpops = psii;
psi_rad_gpops = psii.*pi/180;
gam_gpops = gam;
gam_rad_gpops = gam.*pi/180;
time_gpops = time;

lamX_gpops = lamx;
lamY_gpops = lamy;
lamZ_gpops = lamz;
lamV_gpops = lamv;
lamPSII_gpops = lampsii;
lamGAM_gpops = lamgam;

bank_gpops = bank;
bank_rad_gpops = bank_gpops.*pi/180;
aoa_gpops = aoa;
aoa_rad_gpops = aoa_gpops.*pi/180;
T_gpops = T;
Noise_gpops = Noise;

D_gpops = C1.*v_gpops.^2+C2./(v_gpops.^2);
Ft_gpops = T_gpops.*cos(aoa_rad_gpops) - D_gpops; % force along velocity vector
Fn_gpops = T_gpops.*sin(aoa_rad_gpops) + L; % force perpendicular to velocity vector
 
% Hamiltonian calculation
ham_gpops1 = lamX_gpops.*v_gpops.*cos(gam_rad_gpops).*cos(psi_rad_gpops); 
ham_gpops2 = lamY_gpops.*v_gpops.*cos(gam_rad_gpops).*sin(psi_rad_gpops);
ham_gpops3 = lamZ_gpops.*v_gpops.*sin(gam_rad_gpops);
ham_gpops4 = lamV_gpops.*(Ft_gpops./mass - g.*sin(gam_rad_gpops));
ham_gpops5 = lamPSII_gpops.*(Fn_gpops.*sin(bank_rad_gpops)./(mass.*cos(gam_rad_gpops).*v_gpops));
ham_gpops6 = lamGAM_gpops.*(Fn_gpops.*cos(bank_rad_gpops)./(mass.*v_gpops) - g.*cos(gam_rad_gpops)./v_gpops);
ham_gpops7 = T_gpops.^5.2.*cos(gam_rad_gpops)./(v_gpops.*(z_gpops+50).^2.5);
ham_gpops = ham_gpops1 + ham_gpops2 + ham_gpops3 + ham_gpops4 + ham_gpops5 + ham_gpops6 + ham_gpops7;
cd ..

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%
% Noise and 3DOF trajectory plot
figure(1)
subplot(1,3,1:2)
h1 = plot3(down_bvp4c,cross_bvp4c,alt_bvp4c,'b-*','markersize', 3, 'linewidth', 2);
title('3D Trajectory', 'fontSize', 14 , 'fontWeight' , 'bold')
xlabel('Downrange [km]', 'fontSize', 12 , 'fontWeight' , 'bold')
ylabel('Crossrange [km]', 'fontSize', 12 , 'fontWeight' , 'bold')
zlabel('Altitude [km]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot3(down_gpops,cross_gpops,alt_gpops,'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(1,3,3);
h1 = plot(time_bvp4c,Noise_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Time History Plot for Noise', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Noise [dB]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on;
h2 = plot(time_gpops,Noise_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

% States Time History Plot
figure(2)
subplot(2,3,1)
h1 = plot(time_bvp4c,down_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Time History Plot for Downrange', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Downrange [km]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,down_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,2)
h1 = plot(time_bvp4c,cross_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Time History Plot for Crossrange', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Crossrange [km]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,cross_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,3)
h1 = plot(time_bvp4c,alt_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Time History Plot for Altitude', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Altitude [km]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,alt_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,4)
h1 = plot(time_bvp4c,v_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Time History Plot for Velocity', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Velocity [m/s]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,v_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,5)
h1 = plot(time_bvp4c,psi_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Time History Plot for \psi', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\psi [deg]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,psi_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,6)
h1 = plot(time_bvp4c,gam_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Time History Plot for \gamma', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\gamma [deg]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,gam_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

% Hamiltonian and Control History Plot
figure(3)
subplot(1,3,1)
h1 = plot(time_bvp4c,bank_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Bank Angle History', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Bank Angle [deg]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,bank_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(1,3,2)
h1 = plot(time_bvp4c,aoa_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Angle of Attack History', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Angle of Attack [deg]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,aoa_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(1,3,3)
h1 = plot(time_bvp4c,T_bvp4c, 'b-*','markersize', 3, 'linewidth', 2);
title('Thrust History', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Thrust [N]', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,T_gpops, 'r-o','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

% Costates
figure(4)
subplot(2,3,1)
h1 = plot(time_bvp4c,lamX_bvp4c/1e6, 'b','markersize', 3, 'linewidth', 2);
title('Time History Plot for \lambda_{x}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{x} [dB/m] X 10^{6}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,lamX_gpops/1e6, 'r','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,2)
h1 = plot(time_bvp4c,lamY_bvp4c/1e6, 'b','markersize', 3, 'linewidth', 2);
title('Time History Plot for \lambda_{y}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{y} [dB/m] X 10^{6}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,lamY_gpops/1e6, 'r','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,3)
h1 = plot(time_bvp4c,lamZ_bvp4c/1e10, 'b','markersize', 3, 'linewidth', 2);
title('Time History Plot for \lambda_{z}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{z} [dB/m] X 10^{10}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,lamZ_gpops/1e10, 'r','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,4)
h1 = plot(time_bvp4c,lamV_bvp4c/1e10, 'b','markersize', 3, 'linewidth', 2);
title('Time History Plot for \lambda_{z}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{v} [dB.s/m] X 10^{10}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,lamV_gpops/1e10, 'r','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,5)
h1 = plot(time_bvp4c,lamPSII_bvp4c/1e10, 'b','markersize', 3, 'linewidth', 2);
title('Time History Plot for \lambda_{\psi}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{\psi}[dB/rad] X 10^{10}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,lamPSII_gpops/1e10, 'r','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

subplot(2,3,6)
h1 = plot(time_bvp4c,lamGAM_bvp4c/1e13, 'b','markersize', 3, 'linewidth', 2);
title('Time History Plot for \lambda_{\gamma}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{\gamma} [dB/rad] X 10^{13}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpops,lamGAM_gpops/1e13, 'r','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',12,'FontWeight' , 'bold');
legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 12)

% Hamiltonian
figure(5)
h1 = plot(time_bvp4c,ham_bvp4c, 'b','markersize', 3, 'linewidth', 2);
ylabel('Hamiltonian [dB]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
% h2 = plot(time_gpops,ham_gpops, 'r--','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'BVP4C','GPOPS'},'fontSize', 14)
return