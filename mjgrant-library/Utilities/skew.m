function n = skew(vec)
%#eml
%
% Create a skew symmetric tensor from a vector.
%
% Inputs:
%   vec - a vector of length 3 from which the skew symmetric tensor is created
%
% Outputs:
%   n - resulting skew symmetric tensor
%

% Author: GaTech SSDL / Michael J. Grant on 10/26/08

n = [     0  -vec(3)  vec(2); ...
      vec(3)      0  -vec(1); ...
     -vec(2)  vec(1)      0];

return

