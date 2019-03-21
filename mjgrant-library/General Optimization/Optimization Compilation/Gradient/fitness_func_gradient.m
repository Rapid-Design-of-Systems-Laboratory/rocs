function [out]=fitness_func_gradient(mat_input)

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

if nargin == 0
  out.X = X;
  out.Y = Y;
  out.Z = Z;
  % surf(X,Y,Z);
  % [c,h] = contour(X,Y,Z); colorbar;
else
  out = -Z;
end