function presentation_plot
%
% This function is used to convert current plots to those acceptable for
% presentations.
%

% Author: GaTech SSDL - Michael J. Grant on 11/06/08

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

fontsize_axes = 17-4+4+2; % 15
fontsize_title = 18-4+4+2; % 16
fontsize_axis_labels = 18-4+4+2; % 16

%%%%%%%%%%%%%%%%%%%%
%% Change Numbers %%
%%%%%%%%%%%%%%%%%%%%

h = gca;
set(h,'FontSize',fontsize_axes);

%%%%%%%%%%%%%%%%%%
%% Change Title %%
%%%%%%%%%%%%%%%%%%

h = get(gca,'Title');
set(h,'FontSize',fontsize_title);

%%%%%%%%%%%%%%%%%%%%%%%%
%% Change Axis Labels %%
%%%%%%%%%%%%%%%%%%%%%%%%

h = get(gca,'XLabel');
set(h,'FontSize',fontsize_axis_labels);

h = get(gca,'YLabel');
set(h,'FontSize',fontsize_axis_labels);

h = get(gca,'ZLabel');
set(h,'FontSize',fontsize_axis_labels);

return

