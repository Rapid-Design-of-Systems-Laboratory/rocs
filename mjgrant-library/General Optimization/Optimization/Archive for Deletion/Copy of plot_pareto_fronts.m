function plot_pareto_fronts(mat1,mat2,legend_str,dim1,dim2)
%
% This function is used to compare Pareto fronts.
%

% Author: NASA JSC/DM42 - Michael J. Grant in July 2006


% In Range vs. Alt vs. G directory:
% plot_pareto_fronts('results.mat','results_constrained_2d',{'3D Pareto Surface Projection','2D Pareto Front'},[1 2],[1 2])


% close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data and Build Pareto Front %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load optimization data from mat1
data = load(mat1);
in = data.in;
od = data.od;

% Input logic for mat1
if ~exist('dim1','var') || isempty(dim1)
  dim1 = 1 : 1 : in.num_obj;
end

% Objective information for mat1
obj_label1 = in.obj_label;

% Build Pareto front from mat1
pareto_front1 = build_pareto_front(od.fit(:,:,1:od.iter));

% Load optimization data from mat2 if exists
if exist('mat2','var')
  data = load(mat2);
  in = data.in;
  od = data.od;

  % Input logic regarding mat1
  if ~exist('dim2','var') || isempty(dim2)
    dim2 = 1 : 1 : in.num_obj;
  end

  % Objective information for mat1
  obj_label2 = in.obj_label;

  % Build Pareto front from mat2
  pareto_front2 = build_pareto_front(od.fit(:,:,1:od.iter));
end

%%%%%%%%%%%%%%%%%
%% Error Check %%
%%%%%%%%%%%%%%%%%

if exist('mat2','var')
  % Check dimensions for plotting
  if (length(dim1) > 3) || (length(dim2) > 3)
    error(['dim1 or dim2 is too large. Either error in input or number of ', ...
      'objectives is too large to plot in 3D if no input was specified.']);
  end

  % Check legend cell array
  if exist('legend_str','var') && length(legend) > 2
    error('Legend must be a cell array with length 2.');
  end
end

%%%%%%%%%%%
%% Plots %%
%%%%%%%%%%%

% Plot Pareto front from mat1
if length(dim1) == 2
  plot(pareto_front1(dim1(1),:),pareto_front1(dim1(2),:),'b.');
else
  plot3(pareto_front1(dim1(1),:),pareto_front1(dim1(2),:), ...
    pareto_front1(dim1(3),:),'b.');
end

hold on;
grid on;

if exist('mat2','var')
  % Plot Pareto front from mat2
  if length(dim2) == 2
    plot(pareto_front2(dim2(1),:),pareto_front2(dim2(2),:),'r.');
  else
    plot3(pareto_front2(dim2(1),:),pareto_front2(dim2(2),:), ...
      pareto_front2(dim2(3),:),'r.');
  end

  % Construct legend string
  if ~exist('legend_str','var')
    % Assume legend name
    legend_str = {mat1,mat2};
  end

  % Add legend if legend variable is not empty
  if ~isempty(legend_str)
    legend(legend_str{1},legend_str{2});
  end
end

return

