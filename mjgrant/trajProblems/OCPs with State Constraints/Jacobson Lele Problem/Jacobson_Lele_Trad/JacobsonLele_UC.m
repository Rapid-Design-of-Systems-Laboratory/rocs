% This is the main file for Unconstrained Jaconson Lele Problem
% input : void
% output : void
% Developed by : Kshitij Mall
% Last modified: Mar 8, 2016
close all; clc; clear;
nargin = 0;
if nargin < 1
   solver = 'bvp4c';
end
bvpsolver = fcnchk(solver);

Nt = 100;
tau = linspace(0,1,Nt); % Normalize time variable 

% Guess for time parameters
p(1) = 10; % Time 

% Guess for unknown lambdas/nus
p(2) = 0;
p(3) = 0;
p(4) = 0;
p(5) = 0;
p(6) = 0;
p(7) = 0;
p(8) = 0;
p(9) = 0;
p(10) = 0;
p(11) = 0;

% x1 values
x1 = linspace(0,0,Nt);

% x2 values
x2 = linspace(15/12,-15/12,Nt);

% x3 values
x3 = linspace(-15/12,-15/12,Nt);

% x4 values
x4 = linspace(15/16,-15/16,Nt);

% time values
time = linspace(0,10,Nt);

% lamx1 values
lamx1 = linspace(-0.1,-0.1,Nt);

% lamx2 values
lamx2 = linspace(-0.1,0.1,Nt);

% lamx4 values
lamx3 = linspace(-0.1,0.1,Nt);

% lamx4 values
lamx4 = linspace(-0.1,0.1,Nt);

% lam_time values
lam_time = linspace(-1,1,Nt);

y = [x1;x2;x3;x4;time;lamx1;lamx2;lamx3;lamx4;lam_time];

solinit.x = tau;
solinit.y = y;
solinit.parameters = p; % Pass final time as parameter
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
bvpOptions = bvpset('AbsTol',1e-6,'RelTol',1e-6,'Nmax',100000000,'Stats','on'); % Set options for bvp4c
sol = bvpsolver(@eoms_JacobsonLele_UC,@bcs_JacobsonLele_UC,solinit,bvpOptions);

Z = sol.y;
x1 = Z(1,:); % Altitude, m
x2 = Z(2,:); % Downrange angle, rad
x2 = Z(3,:); % Downrange angle, rad
x2 = Z(4,:); % Downrange angle, rad
time = Z(5,:);
lam_x1 = Z(6,:); % Costate for altitude
lam_x2 = Z(7,:); % Costate for downrange angle
lam_x3 = Z(8,:); % Costate for downrange angle
lam_x4 = Z(9,:); % Costate for downrange angle
lam_t = Z(10,:);
% t = sol.x;
t = time;

u = -lam_x4;

save('results_JLU.mat')

% %%%%%%%%%%
% %% Plot %%
% %%%%%%%%%%
% % Trajectory plot
% figure(1)
% subplot(2,2,1)
% plot(t,x1) 
% title('Displacement history')
% xlabel('Time [s]')
% ylabel('Altitude [m]')
% grid on
% hold on
% 
% % Control history plot
% subplot(2,2,2)
% plot(t,x2)
% title('Velocity history')
% xlabel('Time [s]')
% ylabel('Velocity [m/s]')
% grid on
% hold on
% 
% figure(2)
% subplot(2,2,1)
% plot(t,lam_x1) 
% title('Lamx1 history')
% xlabel('Time [s]')
% ylabel('Lamx1 [m]')
% grid on
% hold on
% 
% % Control history plot
% subplot(2,2,2)
% plot(t,lam_x2)
% title('Lamx2 history')
% xlabel('Time [s]')
% ylabel('Lamx2 [m/s]')
% grid on
% hold on
% 
% figure(3)
% plot(t,u) 
% title('Control history')
% xlabel('Time [s]')
% ylabel('Control, u [m/s^2]')
% grid on
% hold on

% End of file