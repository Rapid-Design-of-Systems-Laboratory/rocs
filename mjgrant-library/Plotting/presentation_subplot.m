function presentation_subplot
%
% This function is used to convert current plots to those acceptable for
% presentations.
%

% Author: GaTech SSDL - Michael J. Grant on 11/06/08

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

fontsize_axes = 17-4+4; % 15
fontsize_title = 18-4+4; % 16
fontsize_axis_labels = 18-4+4; % 16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle Through Subplots %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_set = get(gcf,'Children');

for ctr = 1 : 1 : length(h_set)

%%%%%%%%%%%%%%%%%%%%
%% Change Numbers %%
%%%%%%%%%%%%%%%%%%%%

% h = gca;
set(h_set(ctr),'FontSize',fontsize_axes);

%%%%%%%%%%%%%%%%%%
%% Change Title %%
%%%%%%%%%%%%%%%%%%

h = get(h_set(ctr),'Title');
set(h,'FontSize',fontsize_title);

%%%%%%%%%%%%%%%%%%%%%%%%
%% Change Axis Labels %%
%%%%%%%%%%%%%%%%%%%%%%%%

h = get(h_set(ctr),'XLabel');
set(h,'FontSize',fontsize_axis_labels);

h = get(h_set(ctr),'YLabel');
set(h,'FontSize',fontsize_axis_labels);

end

return

