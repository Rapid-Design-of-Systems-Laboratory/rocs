% Comparison Plot File
clear; close all; clc

% Jacobson Lele Data Trig
load('results_BD_MP_trad.mat')
x1_trig = x1;
x2_trig = x2;
u_trig = u;
time_trig = time;

% Jacobson Lele Data 
load('results_BD_MP.mat')
x1_mp = x1;
x2_mp = x2;
u_mp = u;
time_mp = time;

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

% State x1 Vs Time
figure(1)
subplot(1,2,1)
h1 = plot(time_trig,x1_trig, 'b','markersize', 3, 'linewidth', 4.5);
title('x1 Vs. Time', 'fontSize', 26 , 'fontWeight' , 'bold')
ylabel('x1 [Units]', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_mp,x1_mp, 'r','markersize', 3, 'linewidth', 1.5);
legend([h1 h2],{'Trigonomerized','MPBVP'},'fontSize', 14)
set(gca,'FontSize',24);

% Control Vs Time
subplot(1,2,2)
h3 = plot(time_trig,u_trig, 'b','markersize', 3, 'linewidth', 4.5);
title('Control History', 'fontSize', 26 , 'fontWeight' , 'bold')
ylabel('u [Units]', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h4 = plot(time_mp,u_mp, 'r','markersize', 3, 'linewidth', 1.5);
legend([h3 h4],{'Trigonomerized','MPBVP'},'fontSize', 14)
set(gca,'FontSize',24);

figure(2)
h5 = plot(time_trig,x2_trig, 'b','markersize', 3, 'linewidth', 4.5);
title('Velocity History', 'fontSize', 26 , 'fontWeight' , 'bold')
ylabel('u [Units]', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h6 = plot(time_mp,x2_mp, 'r','markersize', 3, 'linewidth', 1.5);
legend([h5 h6],{'Trigonomerized','MPBVP'},'fontSize', 14)
set(gca,'FontSize',24);
% End of file