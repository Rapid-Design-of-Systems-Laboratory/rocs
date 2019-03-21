function plot_pareto_front(mat,col,dim)
%
% This function is used to compare Pareto fronts.
%

% Author: NASA JSC/DM42 - Michael J. Grant in July 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data and Build Pareto Front %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load optimization data from mat
data = load(mat);
in = data.in;
od = data.od;

% Input logic for mat1
if ~exist('dim','var') || isempty(dim)
  dim = 1 : 1 : in.num_obj;
end

% Build Pareto front from mat
pareto_front = build_pareto_front(od.fit(:,:,1:od.iter));

%%%%%%%%%%%%%%%%%
%% Error Check %%
%%%%%%%%%%%%%%%%%

if length(dim) > 3
  error(['dim1 or dim2 is too large. Either error in input or number of ', ...
      'objectives is too large to plot in 3D if no input was specified.']);
end

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

if ~exist('col','var') || isempty(col)
  col = 'b.';
end

% Plot Pareto front from mat
if length(dim) == 2
  plot(pareto_front(dim(1),:),pareto_front(dim(2),:),col);
else
  plot3(pareto_front(dim(1),:),pareto_front(dim(2),:), ...
    pareto_front(dim(3),:),col);
end

return
