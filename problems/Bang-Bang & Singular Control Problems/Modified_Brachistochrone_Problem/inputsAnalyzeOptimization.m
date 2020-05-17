function [in] = inputsAnalyzeOptimization()

%%%%%%%%%%%%%%%%%
%% Output File %%
%%%%%%%%%%%%%%%%%
% Specify multiple for trade plots

in.outputFile{1} = 'data/results.mat';
ind = 0;

%%%%%%%%%%%%%%%%%%%%%
%% Plotting Inputs %%
%%%%%%%%%%%%%%%%%%%%%

% General plotting items
in.plot.linewidth = 1;
in.plot.markersize = 10;
in.plot.title = true;
in.plot.presentation = true;
in.plot.legend = false;

in.plot.colorSegments.flag = true;
in.plot.colorSegments.colorSet = {'k','r','b'};

% Initial guess
in.plot.initialGuess.flag = false;
in.plot.initialGuess.style = '--';

% Continuation
in.plot.continuation.flag = true;
in.plot.continuation.style = {'k','r','k'};
% in.plot.continuation.style = {'k'};
in.plot.contStartCtr = 0; % 0 = only last solution. Fraction = fraction*numCases
in.plot.contSkip = 1; % How many solutions to skip for plotting, 0 = skip all but last solution (to create envelopes)
in.cont.Index = 2; % Continuation index for plotting

% General movie options
in.movie.fps = 5;
in.movie.quality = 100;
in.movie.compression = 'None';

%%%%%%%%%%%%%%%%%%%
%% Regular Plots %%
%%%%%%%%%%%%%%%%%%%

% Specify regular plots to make (lamR for costate of r, t for time)
% {variable name,bias,scaleFactor,units,plotName}

% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'z','nd','Time'};
% in.figure(ind).plot(1).y{1} = {'x','nd','x'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% 
ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'umax*sin(socontrol)','nd','Control'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% %%%%%%%%%%%%%%
% %% Subplots %%
% %%%%%%%%%%%%%%
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {'lamX','nd','\lambda_1'};
% in.figure(ind).plot(1).legend.location = 'NorthEast';
% in.figure(ind).plot(2).x{1} = {'t','s','Time'};
% in.figure(ind).plot(2).y{1} = {'lamY','nd','\lambda_2'};
% in.figure(ind).plot(2).legend.location = 'NorthWest';
% in.figure(ind).plot(3).x{1} = {'t','s','Time'};
% in.figure(ind).plot(3).y{1} = {'lamALFA','nd','\lambda_{\alpha}'};
% in.figure(ind).plot(3).legend.location = 'NorthWest';
% in.figure(ind).plot(4).x{1} = {'t','s','Time'};
% in.figure(ind).plot(4).y{1} = {'lamX*(V*cos(alfa) + k*epsilon*cos(control)) + lamY*(V*sin(alfa)) + lamALFA*(k*sin(control))','nd','\lambda_{\alpha}'};
% in.figure(ind).plot(4).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
% in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
% in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
% in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];


% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {'x','nd','x'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {'y','nd','y'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
return
