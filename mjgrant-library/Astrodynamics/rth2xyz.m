function [x,y,z] = rth2xyz(r,th,h,omega,theta,inc)
%RTH2XYZ Transform orbital coordinates to Cartesian coordinates
%  [X,Y,Z] = RTH2XYZ(R,TH,H,OMEGA,THETA,INC) transforms cooresponding elements of
%  data stored in orbit-relative coordaintes to Cartesian 
%  coordinates X,Y,Z. The arrays R, TH, and H must be the same size
%  (or any of them can be scalar). The dimensions of X,Y,Z are
%  consistent with the dimensions of R, H.
%
%  Input:
%    R = radius direction component [varying dimension]
%    TH = theta direction component [rad]
%    H = angular momentum direction [varying dimension]
%    OMEGA = longitude of the ascending node [rad]
%    THETA = argument of periapsis + true anomaly [rad]
%    INC = inclination [rad]
%
%  Output:
%    X = x direction component [varying dimension]
%    Y = y direction component [varying dimension]
%    Z = z direction component [varying dimension]
%

%   Michael J. Grant / May 27, 2006


%%%%%%%%%%%%%%%%%%%
%% ASSIGN INPUTS %%
%%%%%%%%%%%%%%%%%%%

co = cos(omega);
so = sin(omega);
ct = cos(theta);
st = sin(theta);
ci = cos(inc);
si = sin(inc);

%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM CONVERSION %%
%%%%%%%%%%%%%%%%%%%%%%%%

x = (co*ct-so*ci*st)*r + (-co*st-so*ci*ct)*th + ( so*si)*h;
y = (so*ct+co*ci*st)*r + (-so*st+co*ci*ct)*th + (-co*si)*h;
z = (si*st)         *r + ( si*ct)         *th + ( ci)   *h;

return

