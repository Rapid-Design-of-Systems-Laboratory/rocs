% This is the main file for Trigonomerized Bryson Denham Problem 
% input : void
% output : void
% Developed by : Kshitij Mall
% Last modified: November 6, 2016

close all; clc; clear;
global x1Max 
x1Max = 0.26;
% scale = 10;
nargin = 0;
if nargin < 1
   solver = 'bvp4c';
end
bvpsolver = fcnchk(solver);

Nt = 1000;
tau = linspace(0,1,Nt); % Normalize time variable for arc 1

% Guess for time parameters
p(1) = 1; % Time for first arc

% Guess for unknown lambdas/nus
p(2) = 0;
p(3) = 0;
p(4) = 0;
p(5) = 0;
p(6) = 0;
p(7) = 0;

% x1 values
x1 = linspace(0,0,Nt);

% x2 values
x2 = linspace(1,-1,Nt);

% time values
time = linspace(0,1,Nt);

% lamx1 values
lamx1 = linspace(-1,1,Nt);

% lamx2 values
lamx2 = linspace(-1,1,Nt);

% lam_time values
lam_time = linspace(-1,1,Nt);

y = [x1;x2;time;lamx1;lamx2;lam_time];

solinit.x = tau;
solinit.y = y;
solinit.parameters = p; % Pass final time as parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
bvpOptions = bvpset('AbsTol',1e-6,'RelTol',1e-6,'Nmax',100000000,'Stats','on'); % Set options for bvp4c
sol = bvpsolver(@eoms_trig,@bcs_trig,solinit,bvpOptions);
% Update
solinit.x = sol.x;
solinit.y = sol.y;
solinit.parameters = sol.parameters; % Pass final time as parameter

Z = sol.y;
x1trig = Z(1,:); % Altitude trig, m
% x1 = x1Max*(sin(x1trig)-epsilon*cos(x1trig));
% x1 = x1Max*2*atan(scale*x1trig)/pi;
x1 = x1Max*sin(x1trig);
x2 = Z(2,:); % Downrange angle, rad
time = Z(3,:);
lam_x1 = Z(4,:); % Costate for altitude
lam_x2 = Z(5,:); % Costate for downrange angle
lam_t = Z(6,:);
% t = sol.x;
t = time;
u = -lam_x2;
%save('results_BD_trig.mat')
%%%%%%%%%%
%% Plot %%
%%%%%%%%%%
% Trajectory plot
figure(1)
subplot(3,2,1)
plot(t,x1,'b*','markersize', 3, 'linewidth', 1) 
xlabel('Time [s]','fontSize', 16 , 'fontWeight' , 'bold')
ylabel('x_{1}','fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
set(gca,'FontSize',18);

subplot(3,2,2)
plot(t,x2,'b*','markersize', 3, 'linewidth', 1)
xlabel('Time [s]','fontSize', 16 , 'fontWeight' , 'bold')
ylabel('x_{2}','fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
set(gca,'FontSize',18);

subplot(3,2,3)
plot(t,lam_x1,'b*','markersize', 3, 'linewidth', 1) 
xlabel('Time [s]','fontSize', 16 , 'fontWeight' , 'bold')
ylabel('\lambda_{x_{1}}','fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
set(gca,'FontSize',18);

subplot(3,2,4)
plot(t,lam_x2,'b*','markersize', 3, 'linewidth', 1)
xlabel('Time [s]','fontSize', 16 , 'fontWeight' , 'bold')
ylabel('\lambda_{x_{2}}','fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
set(gca,'FontSize',18);

subplot(3,2,5)
plot(t,u,'b*','markersize', 3, 'linewidth', 1) 
xlabel('Time [s]','fontSize', 16 , 'fontWeight' , 'bold')
ylabel('u','fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
set(gca,'FontSize',18);
% End of file