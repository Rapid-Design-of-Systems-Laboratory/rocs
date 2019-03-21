% This is the main file for Constrained Bryson Denham Problem
% input : void
% output : void
% Developed by : Kshitij Mall
% Last modified: March 8, 2016
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

% % Guess for time parameters -- This is almost the answer
% p(1) = 0.3748; % Time for first arc
% p(2) = 0.2504; % Time for second arc
% p(3) = 0.3748; % Time for third arc
% 
% % Guess for unknown lambdas/nus
% p(4) = 14.2334;
% p(5) = 5.3354;
% p(6) = 0;
% p(7) = -14.2334;
% p(8) = 5.3354;
% p(9) = 0;
% % Guess for constraint multipliers pi1 and pi2
% p(10) = 28.4668;
% p(11) = 3.5626;

% Guess for time parameters
p(1) = 0.375; % Time for first arc
p(2) = 0.25; % Time for second arc
p(3) = 0.375; % Time for third arc

% Guess for unknown lambdas/nus
p(4) = 0;
p(5) = 0;
p(6) = 0;
p(7) = 0;
p(8) = 0;
p(9) = 0;

% Guess for constraint multipliers pi1 and pi2
p(10) = 0;
p(11) = 0;
p(12) = 0;
% p(12) = 0;
% p(13) = 0;
% p(14) = 0;
% p(15) = 0;

% x1 values
x1_1 = linspace(0,1/8,Nt);
x1_2 = linspace(1/8,1/8,Nt);
x1_3 = linspace(1/8,0,Nt);

% x2 values
x2_1 = linspace(1,0,Nt);
x2_2 = linspace(0,0,Nt);
x2_3 = linspace(0,-1,Nt);

% time values
time_1 = linspace(0,0.375,Nt);
time_2 = linspace(0,0.250,Nt);
time_3 = linspace(0,0.375,Nt);

% lamx1 values
lamx1_1 = linspace(-1,1,Nt);
lamx1_2 = linspace(-1,1,Nt);
lamx1_3 = linspace(-1,1,Nt);

% lamx2 values
lamx2_1 = linspace(-1,1,Nt);
lamx2_2 = linspace(-1,1,Nt);
lamx2_3 = linspace(-1,1,Nt);

% lam_time values
lam_time_1 = linspace(-1,1,Nt);
lam_time_2 = linspace(-1,1,Nt);
lam_time_3 = linspace(-1,1,Nt);

y1 = [x1_1;x2_1;time_1;lamx1_1;lamx2_1;lam_time_1];
y2 = [x1_2;x2_2;time_2;lamx1_2;lamx2_2;lam_time_2];
y3 = [x1_3;x2_3;time_3;lamx1_3;lamx2_3;lam_time_3];

solinit.x = tau;
solinit.y = [y1 y2 y3];
solinit.parameters = p; % Pass final time as parameter
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
bvpOptions = bvpset('AbsTol',1e-6,'RelTol',1e-6,'Nmax',100000000,'Stats','on'); % Set options for bvp4c
sol = bvpsolver(@eoms_trad,@bcs_trad,solinit,bvpOptions);

Z = sol.y;
x1 = Z(1,:); % Altitude, m
x2 = Z(2,:); % Downrange angle, rad
time = Z(3,:);
lam_x1 = Z(4,:); % Costate for altitude
lam_x2 = Z(5,:); % Costate for downrange angle
lam_t = Z(6,:);
% t = sol.x;
t = time;

for count = 1:1:length(sol.x)
    if count<199
        u(count) = -lam_x2(count);
    elseif count>199 && count<300
        u(count) = 0;
    else
        u(count) = -lam_x2(count);
    end
    
end
ham = lam_x1.*x2 + 0.5*u.^2 + lam_t + lam_x2.*u;
save('results_BD_MP_trad.mat')
%%%%%%%%%%
%% Plot %%
%%%%%%%%%%
% Trajectory plot
figure(1)
subplot(2,2,1)
plot(t,x1) 
title('Displacement history')
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on
hold on

% Control history plot
subplot(2,2,2)
plot(t,x2)
title('Velocity history')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on
hold on

figure(2)
subplot(2,2,1)
plot(t,lam_x1) 
title('Lamx1 history')
xlabel('Time [s]')
ylabel('Lamx1 [m]')
grid on
hold on

% Control history plot
subplot(2,2,2)
plot(t,lam_x2)
title('Lamx2 history')
xlabel('Time [s]')
ylabel('Lamx2 [m/s]')
grid on
hold on

figure(3)
plot(t,u) 
title('Control history')
xlabel('Time [s]')
ylabel('Control, u [m/s^2]')
grid on
hold on

% End of file