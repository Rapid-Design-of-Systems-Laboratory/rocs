% Comparison Plot File
clear; close all; clc

% Jacobson Lele Data 
load('results_BD_MP_Trad.mat')
x1_mp = x1;
x2_mp = x2;
u_mp = u;
lam1_mp = lam_x1;
lam2_mp = lam_x2;
ham_mp = ham; 
time_mp = time;

% Jacobson Lele Data 
load('results_BD_MP_Trad.mat')
x1_mp_trad = x1;
x2_mp_trad = x2;
u_mp_trad = u;
lam1_mp_trad = lam_x1;
lam2_mp_trad = lam_x2;
ham_mp_trad = ham; 
time_mp_trad = time;

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

% State x1 Vs Time
figure(1)
subplot(3,2,1)
h1 = plot(time_mp,x1_mp, 'bx','markersize', 3, 'linewidth', 2);
% title('x_1 Vs. Time', 'fontSize', 20 , 'fontWeight' , 'bold')
ylabel('x_1', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_mp_trad,x1_mp_trad, 'r--','markersize', 3, 'linewidth', 1);
legend([h1 h2],{'Enhanced','Traditional'},'fontSize', 14)
set(gca,'FontSize',18);

subplot(3,2,2)
h1 = plot(time_mp,x2_mp, 'bx','markersize', 3, 'linewidth', 2);
% title('x_2 Vs. Time', 'fontSize', 20 , 'fontWeight' , 'bold')
ylabel('x_2', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_mp_trad,x2_mp_trad, 'r--','markersize', 3, 'linewidth', 1);
legend([h1 h2],{'Enhanced','Traditional'},'fontSize', 14)
set(gca,'FontSize',18);

% lam1 Vs Time
subplot(3,2,3)
h1 = plot(time_mp,lam1_mp, 'bx','markersize', 3, 'linewidth', 2);
% title('\lambda_1 Vs. Time', 'fontSize', 20 , 'fontWeight' , 'bold')
ylabel('\lambda_{x_1}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_mp_trad,lam1_mp_trad, 'r--','markersize', 3, 'linewidth', 1);
legend([h1 h2],{'Enhanced','Traditional'},'fontSize', 14)
set(gca,'FontSize',18);

% lam2 Vs Time
subplot(3,2,4)
h1 = plot(time_mp,lam2_mp, 'bx','markersize', 3, 'linewidth', 2);
% title('\lambda_2 Vs. Time', 'fontSize', 20 , 'fontWeight' , 'bold')
ylabel('\lambda_{x_2}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_mp_trad,lam2_mp_trad, 'r--','markersize', 3, 'linewidth', 1);
legend([h1 h2],{'Enhanced','Traditional'},'fontSize', 14)
set(gca,'FontSize',18);

% Control Vs Time
subplot(3,2,5)
h1 = plot(time_mp,u_mp, 'bx','markersize', 3, 'linewidth', 2);
% title('Control History', 'fontSize', 20 , 'fontWeight' , 'bold')
ylabel('u', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_mp_trad,u_mp_trad, 'r--','markersize', 3, 'linewidth', 1);
legend([h1 h2],{'Enhanced','Traditional'},'fontSize', 14)
set(gca,'FontSize',18);

% Control Vs Time
subplot(3,2,6)
h1 = plot(time_mp,ham_mp*1e5, 'bx','markersize', 3, 'linewidth', 2);
% title('Hamiltonian History', 'fontSize', 20 , 'fontWeight' , 'bold')
ylabel('Hamiltonian X 10^{-5}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_mp_trad,ham_mp_trad*1e5, 'r--','markersize', 3, 'linewidth', 1);
legend([h1 h2],{'Enhanced','Traditional'},'fontSize', 14)
set(gca,'FontSize',18);
% End of file