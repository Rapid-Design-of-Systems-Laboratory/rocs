% This is the equations of motion file for Constrained Bryson Denham Problem
% Developed by : Kshitij Mall
% Last modified: 6 Nov, 2016
function dydx = eoms_trig(x,y,p)
global x1Max
% Derivative function 
dydx = zeros(6,1);

% Define the time for different arcs
tSet = p(1);
if y(1) == pi/2

% Based on the control, obtain the time derivatives of states and costates
dydx(1,1) = 0; % y(2);
dydx(2,1) = 0;
dydx(3,1) = 1;
dydx(4,1) = 0; % 0;
dydx(5,1) = 0; % -y(4);
dydx(6,1) = 0; 

else
% Define optimal control based on arcs  
u = -y(5); % For unconstrained arc1

% Based on the control, obtain the time derivatives of states and costates
dydx(1,1) = y(2)/(x1Max*(cos(y(1)))); % y(2);
dydx(2,1) = u;
dydx(3,1) = 1;
dydx(4,1) = -y(4)*sin(y(1))*y(2)/(x1Max*cos(y(1))^2); % 0;
dydx(5,1) = -y(4)/(x1Max*cos(y(1))); % -y(4);
dydx(6,1) = 0;
end

% Based on arcs define the actual total time derivatives
dydx = dydx*tSet;

end

% End of file