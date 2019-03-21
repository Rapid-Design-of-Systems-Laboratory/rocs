% This is the equations of motion file for Constrained Bryson Denham Problem
% Developed by : Kshitij Mall
% Last modified: 18 Mar, 2016
function dydx = eoms_jumps(x,y,region,p)

% Derivative function 
dydx = zeros(6,1);
% Define the time for different arcs
tSet = p(1:4);
% Define optimal control based on arcs
switch region
    case 1    
        u = 1; % For unconstrained arc1
    case 2    
        u = -2*y(5); % For constraint arc2
    case 3
        u = -1; % For unconstrained arc3
    case 4
        u = -2*y(5); % For unconstrained arc3
    otherwise
        error('MATLAB:fourbvp:BadRegionIndex','Incorrect region index: %d',region);
end

% Based on the control, obtain the time derivatives of states and costates
dydx(1,1) = y(2);
dydx(2,1) = -y(1) + y(2)*(1.4 - 0.14*y(2)^2) + 4*u;
dydx(3,1) = 1;
dydx(4,1) = -2*y(1) + y(5);
dydx(5,1) = -y(4) - 1.4*y(5) + 0.42*y(2)^2*y(5);
dydx(6,1) = 0;

% Based on arcs define the actual total time derivatives
dydx = dydx*tSet(region);

end

% End of file