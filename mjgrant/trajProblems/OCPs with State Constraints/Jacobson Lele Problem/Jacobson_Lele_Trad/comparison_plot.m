% Comparison Plot File
clear; close all; clc

% Jacobson Lele Data
load('results_JLU.mat')
x1_u = x1;
u_u = u;
time_u = time;

% Bryson Denham Data
load('results_JL.mat')
ind_const_start = 287;
ind_const_end = 584;

% Time breakdown
time_c1 = time(1:ind_const_start-1);
time_c2 = time(ind_const_start:ind_const_end);
time_c3 = time(ind_const_end+1:end);

% x1 breakdown
x1_c1 = x1(1:ind_const_start-1);
x1_c2 = x1(ind_const_start:ind_const_end);
x1_c3 = x1(ind_const_end+1:end); 

% Control breakdown
u_c1 = u(1:ind_const_start-1);
u_c2 = u(ind_const_start:ind_const_end);
u_c3 = u(ind_const_end+1:end);

x1_c = x1;
u_c = u;
time_c = time;
%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

% State x1 Vs Time
figure(1)
subplot(1,2,1)
h1 = plot(time_u,x1_u, 'r*','markersize', 1, 'linewidth', 4.5);
% title('x1 Vs. Time', 'fontSize', 26 , 'fontWeight' , 'bold')
ylabel('x_{1}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_c,x1_c, 'b--','markersize', 1, 'linewidth', 1.5);
% h2 = plot(time_c1,x1_c1, 'b','markersize', 3, 'linewidth', 1.5);
% h3 = plot(time_c2,x1_c2, 'r','markersize', 3, 'linewidth', 1.5);
% plot(time_c3,x1_c3, 'b','markersize', 3, 'linewidth', 1.5)
% legend([h1 h2 h3],{'Unconstrained Trajectory','Bryson Unconstrained Arc', 'Bryson Constraint Arc'},'fontSize', 14)
legend([h1 h2],{'Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',18);

% Control Vs Time
subplot(1,2,2)
h3 = plot(time_u,u_u, 'r*','markersize', 1, 'linewidth', 4.5);
% title('Control History', 'fontSize', 26 , 'fontWeight' , 'bold')
ylabel('u', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h4 = plot(time_c,u_c, 'b--','markersize', 1, 'linewidth', 1.5);
% h5 = plot(time_c1,u_c1, 'b','markersize', 3, 'linewidth', 1.5);
% h6 = plot(time_c2,u_c2, 'r','markersize', 3, 'linewidth', 1.5);
% plot(time_c3,u_c3, 'b','markersize', 3, 'linewidth', 1.5)
% legend([h4 h5 h6],{'Unconstrained Trajectory','Bryson Unconstrained Arc', 'Bryson Constraint Arc'},'fontSize', 14)
legend([h3 h4],{'Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',18);

% End of file