function [in] = inputsAnalyzeOptimization()

%%%%%%%%%%%%%%%%%
%% Output File %%
%%%%%%%%%%%%%%%%%
% Specify multiple for trade plots

in.outputFile{1} = 'data/del_clh3.mat';
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
in.plot.continuation.style = {'k'};
in.plot.contStartCtr = 1; % 0 = only last solution. Fraction = fraction*numCases
in.plot.contSkip = 0; % How many solutions to skip for plotting, 0 = skip all but last solution (to create envelopes)
in.cont.Index = 3; % Continuation index for plotting

%%%%%%%%%%%%%%%%%%%
%% Regular Plots %%
%%%%%%%%%%%%%%%%%%%

% Specify regular plots to make (lamR for costate of r, t for time)
% {variable name,bias,scaleFactor,units,plotName}
ind = ind+1;
in.figure(ind).plot(1).x{1} = {'v','m/s','Velocity'};
in.figure(ind).plot(1).y{1} = {'h','m','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'thetta*6378','km','Downrange'};
% in.figure(ind).plot(1).y{1} = {'(r-6378000)/1000','km','Altitude'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'gam*180/pi','deg','Flight-Path Angle'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {'CL','deg','Angle of Attack'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'phi*180/pi','deg','Latitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'psii*180/pi','deg','Heading'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'h','m','Altitude'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];


% rho = '(rho0*exp(-(h)/H))'; % Exponential Atmospheric Density [kg/m^3] 
% Cl = '(1.5658*alfa + -0.0000)'; 
% Cd = '(1.6537*alfa^2 + -0.0000*alfa + 0.0612)';
% D = ['(1/2*',rho,'*v^2*',Cd,'*Aref)']; % Drag Force [N] 
% L = ['(1/2*',rho,'*v^2*',Cl,'*Aref)']; % Lift Force [N] 

%%%%%%%%%%%%%%
%% Subplots %%
%%%%%%%%%%%%%%

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'lamH','km/s^2','\lambda_H'};
in.figure(ind).plot(1).legend.location = 'NorthEast';
in.figure(ind).plot(2).x{1} = {'t','s','Time'};
in.figure(ind).plot(2).y{1} = {'lamPHI','km/s','\lambda_v'};
in.figure(ind).plot(2).legend.location = 'SouthEast';
in.figure(ind).plot(3).x{1} = {'t','s','Time'};
in.figure(ind).plot(3).y{1} = {'lamV','km^2/s^2/rad','\lambda_\gamma'};
in.figure(ind).plot(3).legend.location = 'NorthWest';
in.figure(ind).plot(4).x{1} = {'t','s','Time'};
in.figure(ind).plot(4).y{1} = {'lamGAM','km^2/s^2/rad','\lambda_\phi'};
in.figure(ind).plot(4).legend.location = 'NorthEast';
in.figure(ind).plot(5).x{1} = {'t','s','Time'};
in.figure(ind).plot(5).y{1} = {'lamPSII','km^2/s^2/rad','\lambda_\psi'};
in.figure(ind).plot(5).legend.location = 'NorthEast';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'(pi*(1+sin(bankt))/4)*180/pi','deg','Bank Angle'};
in.figure(ind).plot(1).legend.location = 'SouthWest';
in.figure(ind).plot(2).x{1} = {'t','s','Time'};
in.figure(ind).plot(2).y{1} = {'(pi*SA*lamPSII*v*cos((pi*(sin(bankt) + 1))/4)*exp(-h/H)*cos(bankt)*(0.5*CLLB + 0.5*w*(delCLH - (40.4855*(h/hs - 1)^2 + 57.83333*(h/hs - 1)^3 + (13.62131*h)/hs - 10.470043)*(v/(b*h) - (b^2*h^2)/v^2 + (3*b*h)/v - 3) + ((b^2*h^2)/v^2 - (2*b*h)/v + 1)*(10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3 + (2.337815*h)/hs - 1.525574) + ((b^2*h^2)/v^2 - (b*h)/v)*(0.864369*(h/hs - 1)^2 + 12.1*(h/hs - 1)^3 - (2.73417*h)/hs + 3.406847) + (69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3 + (19.0734*h)/hs - 16.705305)*(v^2/(b^2*h^2) - (4*v)/(b*h) + (b^2*h^2)/v^2 - (4*b*h)/v + 6) + (b^2*h^2*(1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3 + (0.834519*h)/hs - 0.723802))/v^2) - 0.5*CLUB*(w - 1) - 0.5*sin(CLw)*(CLLB - w*(delCLH - (40.4855*(h/hs - 1)^2 + 57.83333*(h/hs - 1)^3 + (13.62131*h)/hs - 10.470043)*(v/(b*h) - (b^2*h^2)/v^2 + (3*b*h)/v - 3) + ((b^2*h^2)/v^2 - (2*b*h)/v + 1)*(10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3 + (2.337815*h)/hs - 1.525574) + ((b^2*h^2)/v^2 - (b*h)/v)*(0.864369*(h/hs - 1)^2 + 12.1*(h/hs - 1)^3 - (2.73417*h)/hs + 3.406847) + (69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3 + (19.0734*h)/hs - 16.705305)*(v^2/(b^2*h^2) - (4*v)/(b*h) + (b^2*h^2)/v^2 - (4*b*h)/v + 6) + (b^2*h^2*(1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3 + (0.834519*h)/hs - 0.723802))/v^2) + CLUB*(w - 1))))/(4*cos(gam)) - (pi*SA*lamGAM*v*sin((pi*(sin(bankt) + 1))/4)*exp(-h/H)*cos(bankt)*(0.5*CLLB + 0.5*w*(delCLH - (40.4855*(h/hs - 1)^2 + 57.83333*(h/hs - 1)^3 + (13.62131*h)/hs - 10.470043)*(v/(b*h) - (b^2*h^2)/v^2 + (3*b*h)/v - 3) + ((b^2*h^2)/v^2 - (2*b*h)/v + 1)*(10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3 + (2.337815*h)/hs - 1.525574) + ((b^2*h^2)/v^2 - (b*h)/v)*(0.864369*(h/hs - 1)^2 + 12.1*(h/hs - 1)^3 - (2.73417*h)/hs + 3.406847) + (69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3 + (19.0734*h)/hs - 16.705305)*(v^2/(b^2*h^2) - (4*v)/(b*h) + (b^2*h^2)/v^2 - (4*b*h)/v + 6) + (b^2*h^2*(1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3 + (0.834519*h)/hs - 0.723802))/v^2) - 0.5*CLUB*(w - 1) - 0.5*sin(CLw)*(CLLB - w*(delCLH - (40.4855*(h/hs - 1)^2 + 57.83333*(h/hs - 1)^3 + (13.62131*h)/hs - 10.470043)*(v/(b*h) - (b^2*h^2)/v^2 + (3*b*h)/v - 3) + ((b^2*h^2)/v^2 - (2*b*h)/v + 1)*(10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3 + (2.337815*h)/hs - 1.525574) + ((b^2*h^2)/v^2 - (b*h)/v)*(0.864369*(h/hs - 1)^2 + 12.1*(h/hs - 1)^3 - (2.73417*h)/hs + 3.406847) + (69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3 + (19.0734*h)/hs - 16.705305)*(v^2/(b^2*h^2) - (4*v)/(b*h) + (b^2*h^2)/v^2 - (4*b*h)/v + 6) + (b^2*h^2*(1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3 + (0.834519*h)/hs - 0.723802))/v^2) + CLUB*(w - 1))))/4','nd','Hamiltonian'};
in.figure(ind).plot(2).legend.location = 'NorthEast';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];

return