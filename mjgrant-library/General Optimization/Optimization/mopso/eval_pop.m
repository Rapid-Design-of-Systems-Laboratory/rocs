function [obj,constr] = eval_pop(pos)
%
% [obj] = eval_pop(pos)
%
% This function is used to return all of the objective values and constraint
% values of the population. A simulation could be performed in this function, or
% a continuous function can be tested. If the optimization problem does not have
% any constraints, then constr does not need to be assigned or returned by this
% function.
%
% Input:
%   pos - Position matrix of all individuals in the population. Each individual
%           is represented by a column in pos. Each design variable corresponds
%           to a row of pos.
%
% Output:
%   obj    - Matrix containing the objective values for each individual. The
%            objectives for one individual is formatted in the respective 
%            column.
%   constr - Matrix containing the constraint values for each individual. The
%            constraint values for one individual is formatted in the
%            respective column.
%

% Author: Michael J. Grant - NASA/JSC/DM42 in Feb. 2006


%%%%%%%%%%%%%%%%%%
%% Assign input %%
%%%%%%%%%%%%%%%%%%

X = pos(1,:);
Y = pos(2,:);

%%%%%%%%%%%%%%%%%%%%%%
%% EVALUATE FITNESS %%
%%%%%%%%%%%%%%%%%%%%%%

% % MOP5 (-30 <= x,y <= 30)
% % Journal of Computer Science and Technology Vol. 5, No. 4 (December 2005)
%   Z1 = 0.5*(X.^2 + Y.^2) + sin(X.^2 + Y.^2);
%   Z2 = (3*X - 2*Y + 4).^2/8 + (X - Y + 1).^2/27 + 15;
%   Z3 = 1./(X.^2 + Y.^2 + 1) - 1.1*exp(-X.^2 - Y.^2);

% % MOP6 (0 <= x,y <= 1)
% % Journal of Computer Science and Technology Vol. 5, No. 4 (December 2005)
%   Z1 = X;
%   Z2 = (1 + 10*Y).*(1 - (X./(1 + 10*Y)).^2 - X./(1 + 10*Y).*sin(8*pi*X));

% % MOPC1 (0 <= x <= 5, 0 <= y <= 3)
% % Journal of Computer Science and Technology Vol. 5, No. 4 (December 2005)
%   Z1 = 4*X.^2 + 4*Y.^2;
%   Z2 = (X - 5).^2 + (Y - 5).^2;
% 
%   % Assign constraint values
%   constr(1,:) = (X - 5).^2 + Y.^2 - 25;
%   constr(2,:) = -(X - 8).^2 - (Y + 3).^2 + 7.7;  

% %%%%%%%%%%%%%%%%%%%%%
% %% Test Function A %%
% %%%%%%%%%%%%%%%%%%%%%
% Z1 = -X.^2 - Y.^2;
% Z2 = X + Y + (X+Y).*sin(sqrt(X.^2+Y.^2));

%%%%%%%%%%%%%%%%%%%%%
%% Test Function B %%
%%%%%%%%%%%%%%%%%%%%%
% -5 <= x,y <= 5
Z1 = X + Y;
Z2 = X.^2 + Y.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign Objective Values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Z3','var')
  
  obj = [Z1 ; Z2];
  
elseif in.num_obj == 3
  
  obj = [Z1 ; Z2 ; Z3];
  
end
  
return

