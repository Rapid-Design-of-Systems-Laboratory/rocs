% This is the equations of motion file for Unonstrained Jacobson Lele Problem
% Developed by : Kshitij Mall
% Last modified: 6 Nov, 2016
function dydx = eoms_JacobsonLele_uc(x,y,p)

% Derivative function 
dydx = zeros(10,1);
% Define the time for different arcs
tSet = p(1);
u = -y(9);

% Based on the control, obtain the time derivatives of states and costates
% dydx(1,1) = y(2)/(x1Max*(cos(y(1))+epsilon*sin(y(1))));
dydx(1,1) = y(2);
dydx(2,1) = y(3);
dydx(3,1) = y(4);
dydx(4,1) = u;
dydx(5,1) = 1;
dydx(6,1) = 0;
dydx(7,1) = -y(6);
dydx(8,1) = -y(7);
dydx(9,1) = -y(8);
dydx(10,1) = 0;

% Based on arcs define the actual total time derivatives
dydx = dydx*tSet;

end
% End of file