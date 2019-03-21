function comparison()
close all
cd ./data
load('results_pop1.mat')
data1 = out.setCONT(end).CONT(end).sol;

load('results_pop_xy')
data2 = out.setCONT(end).CONT(end).sol;

load('results_ode')
time_ode = time;
noise_ode = noise_val;
x_ode = x/1000;
y_ode = y/1000;
z_ode = z/1000;
cd ..

% Pop = 1
down_pop1 = data1.y(1,:)/1000;
cross_pop1 = data1.y(2,:)/1000;
alt_pop1 = data1.y(3,:)/1000;
v_pop1 = data1.y(4,:);
psi_pop1 = data1.y(5,:)*180/pi;
gam_pop1 = data1.y(6,:)*180/pi;
tau_pop1 = data1.x;
time_pop1 = data1.parameters(1)*tau_pop1;
lamX_pop1 = data1.y(7,:);
lamY_pop1 = data1.y(8,:);
lamZ_pop1 = data1.y(9,:);
lamV_pop1 = data1.y(10,:);
lamPSII_pop1 = data1.y(11,:);
lamGAM_pop1 = data1.y(12,:);
bank_pop1 = 60*sin(data1.control(1,:));
aoa_pop1 = 15*sin(data1.control(2,:));
T_pop1 = 1560*sin(data1.control(3,:))+1860;
Noise_pop1 = 10*log10(((18.73*T_pop1.^(5.2).*cos(gam_pop1.*pi./180))./(v_pop1.*(alt_pop1.*1000+50).^(2.5))));

% Pop = 5 + y/(x+1)
down_pop_xy = data2.y(1,:)/1000;
cross_pop_xy = data2.y(2,:)/1000;
alt_pop_xy = data2.y(3,:)/1000;
v_pop_xy = data2.y(4,:);
psi_pop_xy = data2.y(5,:)*180/pi;
gam_pop_xy = data2.y(6,:)*180/pi;
tau_pop_xy = data2.x;
time_pop_xy = data2.parameters(1)*tau_pop_xy;
lamX_pop_xy = data2.y(7,:);
lamY_pop_xy = data2.y(8,:);
lamZ_pop_xy = data2.y(9,:);
lamV_pop_xy = data2.y(10,:);
lamPSII_pop_xy = data2.y(11,:);
lamGAM_pop_xy = data2.y(12,:);
bank_pop_xy = 60*sin(data2.control(1,:));
aoa_pop_xy = 15*sin(data2.control(2,:));
T_pop_xy = 1560*sin(data2.control(3,:))+1860;
pop = 4600./(1+data2.y(2,:));
Noise_pop_xy = 10*log10(((18.73*pop.*T_pop_xy.^(5.2).*cos(gam_pop_xy.*pi./180))./(v_pop_xy.*(alt_pop_xy.*1000+50).^(2.5))));

trapz(time_pop1,Noise_pop1)
trapz(time_pop_xy,Noise_pop_xy)
trapz(time_ode,noise_ode)

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%
% Noise and 3DOF trajectory plot
% figure(1)
% subplot(1,3,1:2)
% h1 = plot3(down_pop1,cross_pop1,alt_pop1,'b','markersize', 3, 'linewidth', 2);
% % title('3D Trajectory', 'fontSize', 14 , 'fontWeight' , 'bold')
% xlabel('Downrange [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% ylabel('Crossrange [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% zlabel('Altitude [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot3(down_pop_xy,cross_pop_xy,alt_pop_xy,'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
% subplot(1,3,3)
% h1 = plot(time_pop1,Noise_pop1, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for Noise', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Noise [dB]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,Noise_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
% % States Time History Plot
% figure(2)
% subplot(3,2,1)
% h1 = plot(time_pop1,down_pop1, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for Downrange', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Downrange [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,down_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
% subplot(3,2,2)
% h1 = plot(time_pop1,cross_pop1, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for Crossrange', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Crossrange [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,cross_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
% subplot(3,2,3)
% h1 = plot(time_pop1,alt_pop1, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for Altitude', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Altitude [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,alt_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
% subplot(3,2,4)
% h1 = plot(time_pop1,v_pop1, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for Velocity', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Velocity [m/s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,v_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
% subplot(3,2,5)
% h1 = plot(time_pop1,psi_pop1, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for \psi', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Heading Angle [deg]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,psi_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
% subplot(3,2,6)
% h1 = plot(time_pop1,gam_pop1, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for \gamma', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Flight Path Angle [deg]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,gam_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)
% 
figure(3)
subplot(3,2,1)
h1 = plot(time_pop1,lamX_pop1/1e7, 'b','markersize', 3, 'linewidth', 2);
% title('Time History Plot for \lambda_{x}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{x} X 10^{7} [(dB s)/m]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,lamX_pop_xy/1e7, 'r--','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 12)

subplot(3,2,2)
h1 = plot(time_pop1,lamY_pop1/1e7, 'b','markersize', 3, 'linewidth', 2);
% title('Time History Plot for \lambda_{y}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{y} X 10^{7} [(dB s)/m]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,lamY_pop_xy/1e7, 'r--','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 12)

subplot(3,2,3)
h1 = plot(time_pop1,lamZ_pop1/1e10, 'b','markersize', 3, 'linewidth', 2);
% title('Time History Plot for \lambda_{z}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{z} X 10^{10} [(dB s)/m]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,lamZ_pop_xy/1e10, 'r--','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 12)

subplot(3,2,4)
h1 = plot(time_pop1,lamV_pop1/1e10, 'b','markersize', 3, 'linewidth', 2);
% title('Time History Plot for \lambda_{v}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{v} X 10^{10} [(dB s^{2})/m]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,lamV_pop_xy/1e10, 'r--','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 12)

subplot(3,2,5)
h1 = plot(time_pop1,lamPSII_pop1/1e10, 'b','markersize', 3, 'linewidth', 2);
% title('Time History Plot for \lambda_{\psi}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{\psi} X 10^{10} [(dB s)/rad]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,lamPSII_pop_xy/1e10, 'r--','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 12)

subplot(3,2,6)
h1 = plot(time_pop1,lamGAM_pop1/1e12, 'b','markersize', 3, 'linewidth', 2);
% title('Time History Plot for \lambda_{\gamma}', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{\gamma} X 10^{12} [(dB s)/rad]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,lamGAM_pop_xy/1e12, 'r--','markersize', 3, 'linewidth', 2);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 12)

% Control History Plot
figure(4)
subplot(1,3,1)
h1 = plot(time_pop1,bank_pop1, 'b','markersize', 3, 'linewidth', 2);
% title('Bank Angle History', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Bank Angle [deg]', 'fontSize', 20 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,bank_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)

subplot(1,3,2)
h1 = plot(time_pop1,aoa_pop1, 'b','markersize', 3, 'linewidth', 2);
% title('Angle of Attack History', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Angle of Attack [deg]', 'fontSize', 20 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,aoa_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)

subplot(1,3,3)
h1 = plot(time_pop1,T_pop1/1000, 'b','markersize', 3, 'linewidth', 2);
% title('Thrust History', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('Thrust [kN]', 'fontSize', 20 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 20 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_pop_xy,T_pop_xy/1000, 'r--','markersize', 3, 'linewidth', 1);
set(gca,'FontSize',18,'FontWeight' , 'bold');
legend([h1 h2],{'PDF = 1','PDF = 4600{/}(y+1)'},'fontSize', 14)

% % Noise and 3DOF trajectory plot comparison for optimal trajectory and cda
% figure(5)
% subplot(1,3,1:2)
% h1 = plot3(x_ode,y_ode,z_ode,'b','markersize', 3, 'linewidth', 2);
% % title('3D Trajectory', 'fontSize', 14 , 'fontWeight' , 'bold')
% xlabel('Downrange [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% ylabel('Crossrange [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% zlabel('Altitude [km]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot3(down_pop_xy,cross_pop_xy,alt_pop_xy,'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'CDA','Optimal'},'fontSize', 14)
% 
% subplot(1,3,3)
% h1 = plot(time_ode,noise_ode, 'b','markersize', 3, 'linewidth', 2);
% % title('Time History Plot for Noise', 'fontSize', 14 , 'fontWeight' , 'bold')
% ylabel('Noise [dB]', 'fontSize', 16 , 'fontWeight' , 'bold')
% xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
% grid on
% hold on
% h2 = plot(time_pop_xy,Noise_pop_xy, 'r--','markersize', 3, 'linewidth', 1);
% set(gca,'FontSize',18,'FontWeight' , 'bold');
% legend([h1 h2],{'CDA','Optimal'},'fontSize', 14)
return