% This is the equations of motion file for Constrained Bryson Denham Problem
% Developed by : Kshitij Mall
% Last modified: 17 Mar, 2016
function dydx = eoms_trig(x,y,p)

% Derivative function 
dydx = zeros(6,1);
% Define the time for different arcs
tSet = p(1);

% Define optimal control 
u1 = pi/2;
ham1 = y(4)*y(2) + y(5)*(-y(1) + y(2)*(1.4 - 0.14*y(2)^2) + 4*sin(u1)) + y(6) + sin(u1)^2 + y(1)^2;

u2 = -pi/2;
ham2 = y(4)*y(2) + y(5)*(-y(1) + y(2)*(1.4 - 0.14*y(2)^2) + 4*sin(u2)) + y(6) + sin(u2)^2 + y(1)^2;

if ham1 < ham2
    u = u1;
    hamSave = ham1;
else
    u = u2;
    hamSave = ham2;
end

u3 = asin(-2*y(5));
if imag(u3) == 0
ham3 = y(4)*y(2) + y(5)*(-y(1) + y(2)*(1.4 - 0.14*y(2)^2) + 4*sin(u3)) + y(6) + sin(u3)^2 + y(1)^2;

if ham3<hamSave
    u = u3;
    hamSave = ham3;
end
end

% Based on the control, obtain the time derivatives of states and costates
dydx(1,1) = y(2);
dydx(2,1) = -y(1) + y(2)*(1.4 - 0.14*y(2)^2) + 4*sin(u);
dydx(3,1) = 1;
dydx(4,1) = -2*y(1) + y(5);
dydx(5,1) = -y(4) - 1.4*y(5) + 0.42*y(2)^2*y(5);
dydx(6,1) = 0;

% Based on arcs define the actual total time derivatives
dydx = dydx*tSet;

end

% End of file