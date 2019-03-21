%Driver

clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Target Conditions
alt_target = 0;                 %Altitude Target [m]
vel_target = 0.1;               %Velocity Target [m/s] -- NOTE: Difficult to taget exactly 0 m/s

%Vehicle Parameters
m0 = 300;                  %Vehicle Iniital Mass [kg]
Isp = 218;                      %Specific Impulse [s]
CL = 0.0;                       %Lift Coefficient [-]
CD = 0.0;                       %Drag Coefficient [-]
S = pi/4*10^2;                  %Reference Area [m^2]
Tmax = m0 * 3.718 * 3;          %Maximum Thrust Level [N]

%Planetary Parameters
omega_mag = 0;          %Planetary Rotation Rate [rad/s]
r_planet = 1737.1e3;             %Planetary Radius [m]
mu_planet = 4888376584200;         %Gravitational Parameter [m^3/s^2]

%Initial State
FPA0rel = -60;             %Initial RELATIVE Flight Path Angle [deg]
alt0 = 8660;                    %Initial Altitude [m]
v0rel_mag = 200;           %Initial RELATIVE Velocity [m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Form Initial State
r0 = [ r_planet + alt0; 0 ; 0];
v0rel = v0rel_mag*[sind(FPA0rel);cosd(FPA0rel);0];
v0 = v0rel + cross(omega_mag*[0;0;1],r0);
v0mag = norm(v0)
vtarg = vel_target + omega_mag*(r_planet+alt_target)
deltaV = v0mag - vtarg
X0 = [r0;v0;m0];

%Call Propagator
[m_prop, m_final, T_required, exit_flag] = gravityturn(X0, Isp, CL, CD, S, alt_target,vel_target,Tmax,omega_mag,r_planet,mu_planet)