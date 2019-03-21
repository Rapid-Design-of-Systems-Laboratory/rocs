% This is the main file for Constrained Rayleigh Unconstrained Problem
% input : void
% output : void
% Developed by : Kshitij Mall
% Last modified: March 18, 2016

close all; clc; clear;
nargin = 0;
if nargin < 1
    solver = 'bvp4c';
end
bvpsolver = fcnchk(solver);

Nt = 10;
tau = linspace(0,1,Nt); % Normalize time variable for arc 1

% Guess for time parameters
p(1) = 4.5; % Time for the arc
p(2) = 0;
p(3) = 0;

% Guess for unknown lambdas/nus
p(4) = 0;
p(5) = 0;
p(6) = 0;
p(7) = 0;

% x1 values
x1 = linspace(-5,-0.5,Nt);

% x2 values
x2 = linspace(-5,-1.5,Nt);

% time values
time = linspace(0,4.5,Nt);

% lamx1 values
lamx1 = linspace(-1,1,Nt);

% lamx2 values
lamx2 = linspace(-1,1,Nt);

% lam_time values
lam_time = linspace(0.1,0.1,Nt);

y = [x1;x2;time;lamx1;lamx2;lam_time];

solinit.x = tau;
solinit.y = y;
solinit.parameters = p; % Pass final time as parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
bvpOptions = bvpset('AbsTol',1e-6,'RelTol',1e-6,'Nmax',100000000,'Stats','on'); % Set options for bvp4c
sol = bvpsolver(@eoms_trig,@bcs_trig,solinit,bvpOptions);
toc

Z = sol.y;
x1 = Z(1,:); % Altitude, m
x2 = Z(2,:); % Downrange angle, rad
time = Z(3,:);
lam_x1 = Z(4,:); % Costate for altitude
lam_x2 = Z(5,:); % Costate for downrange angle
lam_t = Z(6,:);
% t = sol.x;
t = time;
u1 = pi/2;
u2 = -pi/2;

for count = 1:1:length(sol.x)
    % Define optimal control
    ham1 = x1(1,count)*x2(1,count) + lam_x2(1,count)*(-x1(1,count) + x2(1,count)*(1.4 - 0.14*x2(1,count)^2) + 4*sin(u1)) + lam_t(1,count) + sin(u1)^2 + x1(1,count)^2;
    ham2 = x1(1,count)*x2(1,count) + lam_x2(1,count)*(-x1(1,count) + x2(1,count)*(1.4 - 0.14*x2(1,count)^2) + 4*sin(u2)) + lam_t(1,count) + sin(u2)^2 + x1(1,count)^2;
    
    if ham1 < ham2
        u(1,count) = u1;
        hamSave = ham1;
    else
        u(1,count) = u2;
        hamSave = ham2;
    end
    
    u3 = (asin(-2*lam_x2(1,count)));
    if imag(u3) == 0
    ham3 = x1(1,count)*x2(1,count) + lam_x2(1,count)*(-x1(1,count) + x2(1,count)*(1.4 - 0.14*x2(1,count)^2) + 4*sin(u3)) + lam_t(1,count) + sin(u3)^2 + x1(1,count)^2;
    
    if ham3<hamSave
        u(1,count) = u3;
    end
    end
end
save('rayleigh_trig.mat');
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
plot(t,sin(u))
title('Control history')
xlabel('Time [s]')
ylabel('Control, u [m/s^2]')
grid on
hold on

% End of file