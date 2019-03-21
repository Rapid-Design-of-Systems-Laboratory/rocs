function [out]=fitness_func(mat_input)
%
% fitness_func is used to return the fitness of all of the test solutions. A
%   simulation could be performed in this function or a continuous function can 
%   be tested. Both the genetic algorithm and the particle swarm utilize this
%   function in order to evaluate the fitness of the indivduals
%
% Input:
%   
%   mat_input -- Assumes that each row specifies a dimension that is varied to
%     find the optimal solution. Thus, each column represents the location of
%     the individual.
%
% Output:
%
%   out -- Matrix containing all of fitness values for each individual
%

if nargin == 0
  % Assumes want to return 3-d plot data to look at function
  x = [-20 : .4 : 20];
  y = [-20 : .4 : 20];
  [X,Y] = meshgrid(x,y);
else
  % Assumes evaluating fitness of test solutions
  X = mat_input(1,:);
  Y = mat_input(2,:);
end

Z1 = -X.^2 - Y.^2;
Z2 = 2*X + Y;
Z3 = 100*sin(X) + 200*sin(Y);
Z = Z1 + Z2 + Z3;
% Z = Z.^2;
% Z = -X.^2-Y.^2;
% Z = X;
if nargin == 0
  out.X = X;
  out.Y = Y;
  out.Z = Z;
  figure;
  surf(X,Y,Z);
  xlabel('x'); ylabel('y'); zlabel('z'); title('Fitness Function');
  figure;
  [c,h] = contour(x,y,Z); colorbar;
  xlabel('x'); ylabel('y'); title('Contour Map of Fitness Function');
else
  out = Z;
end