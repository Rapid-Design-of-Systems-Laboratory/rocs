% This is the equations of motion file for Constrained Bryson Denham Problem
% Developed by : Kshitij Mall
% Last modified: 17 Mar, 2016
function dydx = eoms_uc(x,y,p)

% Derivative function 
dydx = zeros(6,1);
% Define the time for different arcs
tSet = p(1);

% Define optimal control 
% u = -y(5); % For unconstrained arc1
u1 = 0;
ham1 = y(4,end)*y(2,end) + y(5,end)*u1^3 + 0.5*u1^2 + y(6,end);

u2 = -1/(3*y(5,end));
ham2 = y(4,end)*y(2,end) + y(5,end)*u2^3 + 0.5*u2^2 + y(6,end);

if ham1<ham2
    u = u1;
else
    u = u2;
end

% Based on the control, obtain the time derivatives of states and costates
dydx(1,1) = y(2);
dydx(2,1) = u^3; % u
dydx(3,1) = 1;
dydx(4,1) = 0;
dydx(5,1) = -y(4);
dydx(6,1) = 0;

% Based on arcs define the actual total time derivatives
dydx = dydx*tSet;

end

% End of file