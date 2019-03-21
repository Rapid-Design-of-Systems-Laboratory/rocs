% This is the equations of motion file for Constrained Bryson Denham Problem
% Developed by : Kshitij Mall
% Last modified: 8 Feb, 2016
function dydx = eoms_jumps(x,y,region,p)

% Derivative function 
dydx = zeros(6,1);
% Define the time for different arcs
tSet = p(1:3);
% Define optimal control based on arcs
switch region
    case 1    
        u = -y(5); % For unconstrained arc1
    case 2    
        u = 0; % For constraint arc2
    case 3
        u = -y(5); % For unconstrained arc3
    otherwise
        error('MATLAB:fourbvp:BadRegionIndex','Incorrect region index: %d',region);
end

% Based on the control, obtain the time derivatives of states and costates
dydx(1,1) = y(2);
dydx(2,1) = u;
dydx(3,1) = 1;
dydx(4,1) = 0;
dydx(5,1) = -y(4);
dydx(6,1) = 0;

% Based on arcs define the actual total time derivatives
dydx = dydx*tSet(region);

end

% End of file