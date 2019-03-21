% Comparison Plot File
clear; close all; clc

% Jacobson Lele Data Trig
load('results_JLtrig.mat')
x1_trig = x1;
x2_trig = x2;
x3_trig = x3;
x4_trig = x4;
lam_x1_trig = lam_x1;
lam_x2_trig = lam_x2;
lam_x3_trig = lam_x3;
lam_x4_trig = lam_x4;
u_trig = u;
time_trig = time;
ham_trig = lam_x1_trig.*x2_trig+ lam_x2_trig.*x3_trig + lam_x3_trig.*x4_trig + lam_x4_trig.*u_trig + 0.5.*u_trig.^2;

% Jacobson Lele Data 
load('results_JLU.mat')
x1_u = x1;
x2_u = x2;
x3_u = x3;
x4_u = x4;
lam_x1_u = lam_x1;
lam_x2_u = lam_x2;
lam_x3_u = lam_x3;
lam_x4_u = lam_x4;
u_u = u;
time_u = time;
ham_u = lam_x1_u.*x2_u+ lam_x2_u.*x3_u + lam_x3_u.*x4_u + lam_x4_u.*u_u + 0.5.*u_u.^2;
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

% x2 breakdown
x2_c1 = x2(1:ind_const_start-1);
x2_c2 = x2(ind_const_start:ind_const_end);
x2_c3 = x2(ind_const_end+1:end); 

% x3 breakdown
x3_c1 = x3(1:ind_const_start-1);
x3_c2 = x3(ind_const_start:ind_const_end);
x3_c3 = x3(ind_const_end+1:end); 

% x4 breakdown
x4_c1 = x4(1:ind_const_start-1);
x4_c2 = x4(ind_const_start:ind_const_end);
x4_c3 = x4(ind_const_end+1:end); 

% lam_x1 breakdown
lam_x1_c1 = lam_x1(1:ind_const_start-1);
lam_x1_c2 = lam_x1(ind_const_start:ind_const_end);
lam_x1_c3 = lam_x1(ind_const_end+1:end); 

% lam_x2 breakdown
lam_x2_c1 = lam_x2(1:ind_const_start-1);
lam_x2_c2 = lam_x2(ind_const_start:ind_const_end);
lam_x2_c3 = lam_x2(ind_const_end+1:end); 

% lam_x3 breakdown
lam_x3_c1 = lam_x3(1:ind_const_start-1);
lam_x3_c2 = lam_x3(ind_const_start:ind_const_end);
lam_x3_c3 = lam_x3(ind_const_end+1:end); 

% lam_x4 breakdown
lam_x4_c1 = lam_x4(1:ind_const_start-1);
lam_x4_c2 = lam_x4(ind_const_start:ind_const_end);
lam_x4_c3 = lam_x4(ind_const_end+1:end);

% Control breakdown
u_c1 = u(1:ind_const_start-1);
u_c2 = u(ind_const_start:ind_const_end);
u_c3 = u(ind_const_end+1:end);

x1_c = x1;
x2_c = x2;
x3_c = x3;
x4_c = x4;
lam_x1_c = lam_x1;
lam_x2_c = lam_x2;
lam_x3_c = lam_x3;
lam_x4_c = lam_x4;

u_c = u;
mu_c = -lam_x4 - u;
time_c = time;
ham_c = lam_x1_c.*x2_c+ lam_x2_c.*x3_c + lam_x3_c.*x4_c + lam_x4_c.*u_c + 0.5.*u_c.^2 + mu_c.*u;
  
%%%%%%%%%%
%% Plot %%
%%%%%%%%%%
% States Vs Time
figure(1)
subplot(2,2,1)
h1 = plot(time_trig,x1_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('x_{1}', 'fontSize', 24, 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_u,x1_u, 'r--','markersize', 3, 'linewidth', 1.5);
h3 = plot(time_c,x1_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h1 h2 h3],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

subplot(2,2,2)
h4 = plot(time_trig,x2_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('x_{2}', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,x2_u, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,x2_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

subplot(2,2,3)
h4 = plot(time_trig,x3_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('x_{3}', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,x3_u, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,x3_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

subplot(2,2,4)
h4 = plot(time_trig,x4_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('x_{4}', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,x4_u, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,x4_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

% Co-states Vs Time
figure(2)
subplot(2,2,1)
h1 = plot(time_trig,lam_x1_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('\lambda_{x_1} X 10^{-8}', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_u,lam_x1_u, 'r--','markersize', 3, 'linewidth', 1.5);
h3 = plot(time_c,lam_x1_c*1e8, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h1 h2 h3],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

subplot(2,2,2)
h4 = plot(time_trig,lam_x2_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('\lambda_{x_2}', 'fontSize', 24, 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,lam_x2_u, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,lam_x2_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

subplot(2,2,3)
h4 = plot(time_trig,lam_x3_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('\lambda_{x_3}', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,lam_x3_u, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,lam_x3_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

subplot(2,2,4)
h4 = plot(time_trig,lam_x4_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('\lambda_{x_4}', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,lam_x4_u, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,lam_x4_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

figure(3)
% Control Vs Time
subplot(1,2,1)
h4 = plot(time_trig,u_trig, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('u', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,u_u, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,u_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);

% Hamiltonain Vs Time
subplot(1,2,2)
h4 = plot(time_trig,ham_trig*1000, 'b*','markersize', 3, 'linewidth', 1.5);
ylabel('Hamiltonian X 10^{-3}', 'fontSize', 24 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 24 , 'fontWeight' , 'bold')
grid on
hold on
h5 = plot(time_u,ham_u*1000, 'r--','markersize', 3, 'linewidth', 1.5);
h6 = plot(time_c,ham_c, 'k--','markersize', 3, 'linewidth', 1.5);
legend([h4 h5 h6],{'Trigonometrization','Unconstrained Trajectory','Traditional Path Constraint'},'fontSize', 14)
set(gca,'FontSize',20);
% End of file