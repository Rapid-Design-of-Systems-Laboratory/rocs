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
in.plot.colorSegments.colorSet = {'k','r','k'};

% Initial guess
in.plot.initialGuess.flag = false;
in.plot.initialGuess.style = '--';

% Continuation
in.plot.continuation.flag = true;
in.plot.continuation.style = {'k'};
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
ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'h','ft','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'v','ft/s','Velocity'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'mass','kg','Mass'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'96.5*sin(control)+96.5','nd','Thrust'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'lamH','ft','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'lamV','ft','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'lamMASS','ft','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'lamH*(v + epsilon*cos(control)) - lamV*(g0 - (0.5*TMax*(sin(control) + 1) - dragk*v^2*exp(-h/H))/mass) - (0.5*TMax*lamMASS*(sin(control) + 1))/c','Unit','Hamiltonian'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

return

