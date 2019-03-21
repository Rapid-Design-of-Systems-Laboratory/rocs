function [v1,v2] = gauss(r1,r2,dt,type,mu,tol,output)
%
%  [v1,v2] = gauss(r1,r2,dt,type,mu,tol,output)
%
%  Solves Gauss' problem using universal variables. The velocity vectors are
%  determined at each position vector given the time of flight and direction of
%  motion.
%
% Input:
%   r1 = 3x1 vector of position 1 [varies]
%   r2 = 3x1 vector of position 2 [varies]
%   dt = time of flight [varies]
%   type = specifies direction of motion
%     'short' - short-way transfer
%     'long' - long-way transfer
%   mu = gravitational parameter [varies]
%   tol = Newtonian iteration tolerance on time of flight [varies]
%   output = determines if output should be displayed
%     'true'  = display output
%     'false' = suppress output
%
% Output:
%   v1 = 3x1 vector of velocity at position 1 [varies]
%   v2 = 3x1 vector of velocity at position 2 [varies]
%

% Author: Michael J. Grant / Oct. 20, 2006

%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

if output
  fprintf('\n');
  fprintf('%5s %15s %15s %15s %15s %15s %15s\n', ...
    'i','z','y','x','tc','dt/dz','z_new');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1 - Determine Transfer Angle and A %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute transfer angle
if strcmp(type,'short')
  % Short-way transfer. Transfer angle < pi.
  dTheta = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
elseif strcmp(type,'long')
  % Long-way transfer. Transfer angle > pi.
  dTheta = 2*pi - acos(dot(r1,r2)/(norm(r1)*norm(r2)));
end

% Compute A (constant)
A = sqrt(norm(r1)*norm(r2))*sin(dTheta)/sqrt(1-cos(dTheta));

% If short-way transfer, compute minimum z value allowed
if strcmp(type,'short')
  z_min = fzero(@(z) find_zmin(z,A,norm(r1),norm(r2)),0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2 - Choose Initial z Value %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not change! This initial guess is valuable for problems in which r1 and r2
% are nearly colinear. The other singularity occurs when the TOF approaches
% zero. This scenario (x = 0) is accounted for later in the code.
z =  dTheta^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Iterative Search of z %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;
while 1

  % Iteration number
  i = i + 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Step 3 - Evaluate S(z) and C(z) %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % The following equations are mathematically identical when taking into
  % consideration imaginary numbers. Separate equations (below) used for 
  % computational purposes.
  if z > 0
    C = 1/z*(1-cos(sqrt(z)));
    S = 1/z^(3/2)*(sqrt(z) - sin(sqrt(z)));
  elseif z == 0
    C = 1/2;
    S = 1/6;
  elseif z < 0
    C = (1 - cosh(sqrt(-z)))/z;
    S = (sinh(sqrt(-z)) - sqrt(-z))/(-z)^(3/2);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Step 4 - Compute Auxilliary Variable y and x %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Auxilliary variable y
  y = norm(r1) + norm(r2) - A/sqrt(C)*(1-z*S);
  
  % Universal variable x
  x = sqrt(y/C);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Step 5 - Compute Time of Flight and Compare %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Compute time of flight for given z
  dt_i = (A*sqrt(y) + (y/C)^(3/2)*S)/sqrt(mu);
  
  % Compare with tolerance
  if abs(dt_i - dt) < tol
    % Solution converged
    break;
  else
    % Solution not converged. Perform Newtonian iteration.
    
    

    % Compute dt/dz (if applicable) and new z
    z_old = z;
    if x == 0
      
      % x = 0 => TOF = 0. Set z to small value above z_min (which corresponds to
      % 0 TOF).
      z = z_min + abs(z_min)*1e-12;
      
    else
      
      % Compute dC/dz
      if z ~= 0
        % z does not equal zero. Use usual C and S derivatives.
        dC = 1/(2*z)*(1 - z*S - 2*C);
        dS = 1/(2*z)*(C - 3*S);
      else
        % z equals zero. Usual C and S derivatives break down. Use Taylor Series.
        dC = 1/24;
        dS = 1/120;
      end
      
      % Compute dt/dz
      dt_dz = 1/sqrt(mu)*(x^3*(dS - 3*S*dC/(2*C)) + A/8*(3*S*sqrt(y)/C + A/x));
      
      % Compute new z
      z = z_old + (dt - dt_i)/dt_dz;
      
      % Assumes vehicle does not orbit more than once
      if z > (2*pi)^2
        z = (2*pi)^2-0.001;
      end
      
    end
    
    % If new z is less than z_min (for short transfers), set to z_min
    if strcmp(type,'short') && (z < z_min)
      z = z_min;
    end
    
    % Output Newtonian iteration info to screen
    if output
      fprintf('%5i %15g %15g %15g %15g %15g %15g\n', ...
        i,z_old,y,x,dt_i,dt_dz,z);
    end
    
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time of Flight Converged %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 6 - Compute f and g %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = 1 - y/norm(r1);
g = A*sqrt(y/mu);
gdot = 1 - y/norm(r2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 7 - Compute v1 and v2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1 = (r2 - f*r1)/g;
v2 = (gdot*r2 - r1)/g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 8 - Perform Checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute energy
energy1 = norm(v1)^2/2 - mu/norm(r1);
energy2 = norm(v2)^2/2 - mu/norm(r2);

% Energy should remain constant at both positions
if abs(energy1 - energy2) > 1e-3
  error('Energy not constant at r1 and r2!');
end

% Compute angular momentum
h1 = norm(cross(r1,v1));
h2 = norm(cross(r2,v2));

% Angular momentum should remain constant at both positions
if abs(h1 - h2) > 1e-3
  error('Angular momentum not constant at r1 and r2!');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = find_zmin(z,A,r1,r2)
%
% The function is used to determine the minimum z value permitted for short-way
% transfers.
%
% Input:
%   z = Universal variable z
%   A = Constant value
%   r1 = Magnitude of position vector 1
%   r2 = Magnitude of position vector 2
%
% Output:
%   f = function value in which a zero value corresponds to the minimum z value
%   allowed for short-way transfers
%

% Author: Michael J. Grant / Oct. 20, 2006

% Compute C and S
if z > 0
  C = 1/z*(1-cos(sqrt(z)));
  S = 1/z^(3/2)*(sqrt(z) - sin(sqrt(z)));
elseif z == 0
  C = 1/2;
  S = 1/6;
elseif z < 0
  C = (1 - cosh(sqrt(-z)))/z;
  S = (sinh(sqrt(-z)) - sqrt(-z))/(-z)^(3/2);
end

% Compute function
f = A/sqrt(C)*(1-z*S) - r1 - r2;

return

