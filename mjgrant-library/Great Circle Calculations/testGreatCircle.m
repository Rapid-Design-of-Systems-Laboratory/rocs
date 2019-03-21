function testGreatCircle
%
% Assumes flight < 1/2 circumference of spherical Earth.
%
% A: initial point
% B: target point
% C: point on arc that is shortest distance from current point
% D: current point

close all; clc;

%% Inputs

% Planetary constants
rMag = (6.3781e6+6.3568e6)/2;

% Initial point
A.lat = 0*pi/180;
A.lon = 0*pi/180;

% Target point
B.lat = 0*pi/180;
B.lon = 50*pi/180;

% Current point
D.lat = 20*pi/180;
D.lon = 30*pi/180;

[rangeGo,rangeFlown,xRange] = greatCircleCalcs(A,B,D,rMag)

return

