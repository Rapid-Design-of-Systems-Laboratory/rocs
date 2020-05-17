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
in.cont.Index = 10; % Continuation index for plotting

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
in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
in.figure(ind).plot(1).y{1} = {'h/1000','km','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'thetta*6378','km','Downrange'};
in.figure(ind).plot(1).y{1} = {'h/1000','km','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'gam*180/pi','deg','Flight-Path Angle'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'40*sin(alfamctrig)','deg','Angle of Attack'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'sqrt((1/2*(rho0*exp(-h/H))*v^2*(Cl1*alfamax*sin(alfamctrig))*Aref)^2 + (1/2*(rho0*exp(-h/H))*v^2*(Cd2*alfamax*sin(alfamctrig)^2 + Cd0)*Aref)^2)/(mass*g)','nd','g-load'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'k*sqrt((rho0*exp(-h/H))/rn)*v^3*1e-4','W/cm^2','Heat Rate'};
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
% in.figure(ind).plot(1).y{1} = {'lamH','s/m','\lambda_R'};
% in.figure(ind).plot(1).legend.location = 'NorthEast';
% in.figure(ind).plot(2).x{1} = {'t','s','Time'};
% in.figure(ind).plot(2).y{1} = {'lamTHETTA','s/rad','\lambda_\theta'};
% in.figure(ind).plot(2).legend.location = 'NorthWest';
% in.figure(ind).plot(3).x{1} = {'t','s','Time'};
% in.figure(ind).plot(3).y{1} = {'lamV','s^2/m','\lambda_v'};
% in.figure(ind).plot(3).legend.location = 'SouthEast';
% in.figure(ind).plot(4).x{1} = {'t','s','Time'};
% in.figure(ind).plot(4).y{1} = {'lamGAM','s/rad','\lambda_\gamma'};
% in.figure(ind).plot(4).legend.location = 'NorthWest';
% 
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
% in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
% in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
% in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];

return

