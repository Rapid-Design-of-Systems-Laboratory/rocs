function run_gravity_turn_edl_notes

close all; clc;

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

g = 1.62; % Gravity (assume constant) [m/s^2]
ge = 9.81; % Gravity acceleration at Earth's surface [m/s^2]
Isp = 218; % Specific impulse [s]

% % EDL Notes Example
% V0 = 200; % Velocity
% gamma0 = -60*pi/180; % Flight path angle
% h0 = 8660; % Altitude
% x0 = 0; % Downrange
% m0 = 300; % Mass

% My case
V0 = 180;
gamma0 = -25*pi/180;
h0 = 11000;
x0 = 0;
m0 = 2000;
Isp = 194;
g = 3.7176;

%%%%%%%%%%%%%%%%%%%%%%
%% Run Gravity Turn %%
%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions
y0 = [V0; gamma0; h0; x0; m0];

% Assign parameters
param.g = g;
param.ge = ge;
param.Isp = Isp;

% Fly until reach just above ground (get singularity of guidance at ground)
options = odeset('AbsTol',1e-10,'RelTol',1e-10, ...
  'Events',@gravity_turn_event_func);
[t,y] = ode45(@gravity_turn_eom,[0 1e5],y0,options,param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = y(:,1); % Velocity
gamma = y(:,2); % Flight path angle
h = y(:,3); % Altitude
x = y(:,4); % Downrange
m = y(:,5); % Mass

R = h./sin(-gamma); % Slant range
T = m.*(g+V.^2./(2*R)); % Required thrust
n = T./(m*g);
n1 = 1 + V.^2./(2*R*g);

%%%%%%%%%%%
%% Plots %%
%%%%%%%%%%%

figure;
plot(t,h);
xlabel('Time [sec]');
ylabel('Altitude [sec]');
title('Altitude vs. Time');
grid on;

figure;
plot(t,gamma*180/pi);
xlabel('Time [sec]');
ylabel('Flight Path Angle [deg]');
title('Flight Path Angle vs. Time');
grid on;

figure;
plot(t,V);
xlabel('Time [sec]');
ylabel('Velocity [m/s]');
title('Velocity vs. Time');
grid on;

figure;
plot(t,T);
xlabel('Time [sec]');
ylabel('Thrust [N]');
title('Thrust vs. Time');
grid on;

figure;
plot(t,n);
hold on;
xlabel('Time [sec]');
ylabel('G-Loading');
title('G-Loading vs. Time');
grid on;

pmr = (m(1)-m(end))/m(1)

return

