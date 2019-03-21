function state = gaplotbestfun(options,state,flag)
%GAPLOTBESTF Plots the best score and the mean score.
%   STATE = PLOTBEST(OPTIONS,STATE,FLAG) plots the best score as well
%   as the mean of the scores.
%
%   Example:
%    Create an options structure that will use PLOTBEST
%    as the plot function
%     options = gaoptimset('PlotFcns',@plotbest);

%   Copyright 2004 The MathWorks, Inc. 
%   $Revision: 1.9 $  $Date: 2004/01/16 16:51:50 $

if(strcmp(flag,'init'))
    set(gca,'xlim',[1,options.Generations]);
    xlabel('Generation','interp','none');
    grid on
    ylabel('Fitness value','interp','none');
end

hold on;
generation = state.Generation;
best = min(state.Score);
plot(generation,best, 'v');
title(['Best: ',num2str(best)],'interp','none')
hold off;
        