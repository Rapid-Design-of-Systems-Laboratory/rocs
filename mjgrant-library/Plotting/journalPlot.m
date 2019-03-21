function journalPlot
%
% This function is used to convert current plots to those acceptable for
% presentations.
%

% Author: GaTech SSDL - Michael J. Grant on 11/06/08

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

numPlots = 2;
mapColor = colormap(gray(numPlots+1));
markerStyles = {'o','s','^'}

% fontsize_axes = 13;
% % fontsize_title = 14;
% fontsize_axis_labels = 14;

%%%%%%%%%%%%%%%%%%%%%%
%% Change Font Size %%
%%%%%%%%%%%%%%%%%%%%%%

presentation_subplot;

%%%%%%%%%%%%%%%%%
%% Change Size %%
%%%%%%%%%%%%%%%%%

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 12 5]); % 3.25 3.25 or 12 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change Color and Line Style %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_set = findall(gcf,'type','line');
plotChangeInd = 0;
markerChangeInd = 0;
parentSet = [];

for ctr = 1 : 1 : length(h_set)
  
  h = h_set(ctr);
  if length(get(h,'XData')) > 3
    % Change color to grayscale
    plotChangeInd = plotChangeInd + 1;
    set(h,'Color',mapColor(plotChangeInd,:));
    get(h)
    
    % Change marker if one used
    if ~strcmp(get(h,'Marker'),'none')
      markerChangeInd = markerChangeInd + 1;
      set(h,'Marker',markerStyles{markerChangeInd});
    end
    
    % Get parent to see if in same axes
    hParent = get(h,'Parent');
    if ~ismember(hParent,parentSet)
      parentSet = [parentSet; hParent];
      plotChangeInd = 0;
      markerChangeInd = 0;
    end
    
  end
  
%   if plotChangeInd == numPlots
%     plotChangeInd = 0;
%   end
%   
%   if markerChangeInd == length(markerStyles)
%     markerChangeInd = 0;
%   end
  
% h = gca;
% h = get(h_set(ctr),'Children')
% h = get(gca,'Children');
% x = get(h,'XData')
% y = get(h,'YData')
  
%   set(h_set(ctr),'color',[1 0 0]);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Figure as EPS %%
%%%%%%%%%%%%%%%%%%%%%%%%

% saveas(gcf,'output','eps');

print -depsc2 -r600 outputPlot

return

