% This is the main file for Constrained Jacobson Lele Problem
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
tau1 = linspace(0,1,Nt); % Normalize time variable for arc 1
tau2 = 1+tau1; % Normalize time variable for arc 2
tau3 = 1+tau2; % Normalize time variable for arc 3
tau = [tau1 tau2 tau3];

% Guess for time parameters
p(1) = 3.5; % Time for first arc
p(2) = 3; % Time for second arc
p(3) = 3.5; % Time for third arc

% Guess for unknown lambdas/nus
p(4) = 0;
p(5) = 0;
p(6) = 0;
p(7) = 0;
p(8) = 0;
p(9) = 0;
p(10) = 0;
p(11) = 0;
p(12) = 0;
p(13) = 0;

% Guess for constraint multipliers pi1, pi2, pi3, pi4
p(14) = 0;
p(15) = 0;
p(16) = 0;
p(17) = 0;

% x1 values
x1_1 = linspace(0,1,Nt);
x1_2 = linspace(1,1,Nt);
x1_3 = linspace(1,0,Nt);

% x2 values
x2_1 = linspace(15/12,0,Nt);
x2_2 = linspace(0,0,Nt);
x2_3 = linspace(0,-15/12,Nt);

% x3 values
x3_1 = linspace(-15/12,0,Nt);
x3_2 = linspace(0,0,Nt);
x3_3 = linspace(0,-15/12,Nt);

% x4 values
x4_1 = linspace(15/16,0,Nt);
x4_2 = linspace(0,0,Nt);
x4_3 = linspace(0,-15/16,Nt);

% time values
time_1 = linspace(0,3.75,Nt);
time_2 = linspace(0,2.50,Nt);
time_3 = linspace(0,3.75,Nt);

% lamx1 values
lamx1_1 = linspace(-0.1,0.1,Nt);
lamx1_2 = linspace(-0.1,0.1,Nt);
lamx1_3 = linspace(-0.1,0.1,Nt);

% lamx2 values
lamx2_1 = linspace(-0.1,0.1,Nt);
lamx2_2 = linspace(0,0,Nt);
lamx2_3 = linspace(-0.1,0.1,Nt);

% lamx3 values
lamx3_1 = linspace(-0.1,0.1,Nt);
lamx3_2 = linspace(0,0,Nt);
lamx3_3 = linspace(-0.1,0.1,Nt);

% lamx4 values
lamx4_1 = linspace(-0.1,0.1,Nt);
lamx4_2 = linspace(0,0,Nt);
lamx4_3 = linspace(-0.1,0.1,Nt);

% lam_time values
lam_time_1 = linspace(-0.1,0.1,Nt);
lam_time_2 = linspace(0,0,Nt);
lam_time_3 = linspace(-0.1,0.1,Nt);

y1 = [x1_1;x2_1;x3_1;x4_1;time_1;lamx1_1;lamx2_1;lamx3_1;lamx4_1;lam_time_1];
y2 = [x1_2;x2_2;x3_2;x4_2;time_2;lamx1_2;lamx2_2;lamx3_2;lamx4_2;lam_time_2];
y3 = [x1_3;x2_3;x3_3;x4_3;time_3;lamx1_3;lamx2_3;lamx3_3;lamx4_3;lam_time_3];

solinit.x = tau;
solinit.y = [y1 y2 y3];
solinit.parameters = p; % Pass final time as parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
bvpOptions = bvpset('AbsTol',1e-6,'RelTol',1e-6,'Nmax',100000000,'Stats','on'); % Set options for bvp4c
sol = bvpsolver(@eoms_JacobsonLele,@bcs_JacobsonLele,solinit,bvpOptions);

Z = sol.y;
x1 = Z(1,:); % Altitude, m
x2 = Z(2,:); % Downrange angle, rad
x3 = Z(3,:); % Downrange angle, rad
x4 = Z(4,:); % Downrange angle, rad
time = Z(5,:);
lam_x1 = Z(6,:); % Costate for altitude
lam_x2 = Z(7,:); % Costate for downrange angle
lam_x3 = Z(8,:); % Costate for downrange angle
lam_x4 = Z(9,:); % Costate for downrange angle
lam_t = Z(10,:);
% t = sol.x;
t = time;

for count = 1:1:length(sol.x)
    if count<287
        u(count) = -lam_x4(count);
    elseif count>=287 && count<=584
        u(count) = 0;
    else
        u(count) = -lam_x4(count);
    end
    
end

% keyboard
save('results_JL.mat')

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
% % figure(3)
% % plot(t,u) 
% % title('Control history')
% % xlabel('Time [s]')
% % ylabel('Control, u [m/s^2]')
% % grid on
% % hold on

% End of file