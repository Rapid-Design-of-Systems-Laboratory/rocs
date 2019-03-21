% greatCircleCalcs.m 
%
% Downrange and crossrange calculations relative to great circle connecting
% initial and target surface locations.
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   A - data structure, rad, initial location
%       A.lat, double(1), rad, initial latitude
%       A.lon, double(1), rad, initial longitude
%   B - data structure, rad, target location
%       B.lat, double(1), rad, target latitude
%       B.lon, double(1), rad, target longitude
%   D - data structure, rad, current location
%       B.lat, double(1), rad, current latitude
%       B.lon, double(1), rad, current longitude
%   rMag - double(1), m, radius of assumed spherical planet
%   
% Outputs:
%   rangeGo - double(1), m, downrange to target
%   rangeFlown - double(1), m, downrange from initial location
%   xRange - double(1), m, crossrange from great circle connecting initial
%     and target locations
%
% Major Revision History:
%   *Created 2012, M. Grant
%

function [rangeGo,rangeFlown,xRange] = greatCircleCalcs(A,B,D,rMag)
%#codegen

%% Convert to Cartesian position vectors (normalized)
Ahat = [cos(A.lat)*cos(A.lon); cos(A.lat)*sin(A.lon); sin(A.lat)];
Bhat = [cos(B.lat)*cos(B.lon); cos(B.lat)*sin(B.lon); sin(B.lat)];
Dhat = [cos(D.lat)*cos(D.lon); cos(D.lat)*sin(D.lon); sin(D.lat)];

%% Compute normal vector of plane containing A and B vectors
nHat = cross(Ahat,Bhat)/norm(cross(Ahat,Bhat));

%% Compute distance from current point to target point
dist2Plane = rMag*dot(Dhat,nHat);

%% Determine Cartesian position of projected vector on plane
Cproj = rMag*Dhat - dist2Plane*nHat;

%% Determine point on great circle that corresponds to this projected vector
Chat = Cproj/norm(Cproj);

%% Determine distance between current point and projected point on circular arc (crossrange angle)
angXrange = acos(dot(Chat,Dhat));
xRange = angXrange*rMag;

%% Determine angle between projected point on circular arc and target (downrange to go)
angRangeGo = acos(dot(Chat,Bhat));
rangeGo = angRangeGo*rMag;

%% Determine angle between projected point on circular arc and initial point (downrange flown)
angRangeFlown = acos(dot(Chat,Ahat));
rangeFlown = angRangeFlown*rMag;

return

