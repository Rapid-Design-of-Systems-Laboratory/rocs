function [ydot]= gravity_turn_eom(t,y,param)

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

V = y(1); % Velocity
gamma = y(2); % Flight path angle
h = y(3); % Altitude
x = y(4); % Downrange
m = y(5); % Mass

g = param.g;
ge = param.ge;
Isp = param.Isp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Guidance Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = h/sin(-gamma); % Slant range
T = m*(g+V^2/(2*R)); % Required thrust

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derivative Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t
V
h
dV_dt = -T/m - g*sin(gamma) - 1/2*0.001*V^2*1.6*pi/4*(4.5^2); % d(velocity)/d(time)
dgamma_dt = -g/V*cos(gamma); % d(flight path angle)/d(time)
dh_dt = V*sin(gamma); % d(altitude)/d(time)
dx_dt = V*cos(gamma); % d(downrange)/d(time)
dm_dt = -T/(ge*Isp); % d(mass)/d(time)

ydot = [dV_dt; dgamma_dt; dh_dt; dx_dt; dm_dt];

return

