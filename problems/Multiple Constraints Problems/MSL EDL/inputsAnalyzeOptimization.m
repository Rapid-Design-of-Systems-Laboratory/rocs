function [in] = inputsAnalyzeOptimization()

%%%%%%%%%%%%%%%%%
%% Output File %%
%%%%%%%%%%%%%%%%%

in.outputFile{1} = 'data/msledl_bvp4c_ac.mat';
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
in.plot.contStartCtr = 0; % 0 = only last solution. Fraction = fraction*numCases
in.plot.contSkip = 1; % How many solutions to skip for plotting, 0 = skip all but last solution (to create envelopes)
in.cont.Index = 6; % Continuation index for plotting

% General movie options
in.movie.fps = 5;
in.movie.quality = 100;
in.movie.compression = 'None';

%%%%%%%%%%%%%%%%%%%
%% Regular Plots %%
%%%%%%%%%%%%%%%%%%%

% % Specify regular plots to make (lamR for costate of r, t for time)
% % {variable name,bias,scaleFactor,units,plotName}
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'(r - a*(1-0.5*(1-cos(2*thetta))*((a-b)/a) - (a/(4*r)-1/16)*(1-cos(4*thetta))*((a-b)/a)^2))/1000','km','Altitude'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
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

% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'acosd(0.683*sin(sigmactrl) - 0.183)','deg','Control'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'t','s','Time'};
in.figure(ind).plot(1).y{1} = {'acosd(0.683*sin(sigmaaactrl) + 0.183)','deg','Control'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% Heat Rate
ind = ind+1;
in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
in.figure(ind).plot(1).y{1} = {'(k*sqrt((rho0*exp(-h/H))/rn)*v^3)/1e4','W/cm^2','Heat Rate'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% G-Load
ind = ind+1;
in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
in.figure(ind).plot(1).y{1} = {'(sqrt((1/2*(rho0*exp(-h/H))*v^2*Cd*LbyD*Aref)^2 + (1/2*(rho0*exp(-h/H))*v^2*Cd*Aref)^2)/(mass*g))','Earth-gs','G-Load'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% Dynamic Pressure
ind = ind+1;
in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
in.figure(ind).plot(1).y{1} = {'(0.5*(rho0*exp(-h/H))*v^2)/1000','kPa','Dynamic Pressure'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% 
ind = ind+1;
in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
in.figure(ind).plot(1).y{1} = {'epsQ/cos((0.5*k*v^3*pi*((rho0*exp(-h/H))/rn)^(1/2))/Qmax) - lamV*((muu*sin(gam))/(h + rm)^2 + (Aref*Cd*rho0*v^2*exp(-h/H))/(2*mass)) + epsq/cos((0.25*rho0*v^2*pi*exp(-h/H))/qmax) + epsg/cos((0.5*pi*((exp(-(2*h)/H)*Aref^2*Cd^2*LbyD^2*rho0^2*v^4)/4 + (exp(-(2*h)/H)*Aref^2*Cd^2*rho0^2*v^4)/4)^(1/2))/(g*gmax*mass)) + lamGAM*((v*cos(gam))/(h + rm) - (muu*cos(gam))/(v*(h + rm)^2) + (Aref*Cd*LbyD*rho0*v*exp(-h/H)*(0.683*sin(sigmaaactrl) + 0.183))/(2*mass)) + epsc*cos(sigmaaactrl) + lamH*v*sin(gam)','Hamiltonian','Hamiltonian'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

ind = ind+1;
in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
in.figure(ind).plot(1).y{1} = {'epsc*cos(sigmaaactrl)','m/s','ErrorControlPart'};
in.figure(ind).plot(1).legend.location = 'NorthWest';
in.figure(ind).movie.make = false;
in.figure(ind).movie.name = 'testMovie';
in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'epsq*sec(0.5*pi*(2*(0.5*(rho0*exp(-h/H))*v^2)-qmax)/qmax)','m/s','qpart'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'epsg*sec(0.5*pi*(2*(sqrt((1/2*(rho0*exp(-h/H))*v^2*Cd*LbyD*Aref)^2 + (1/2*(rho0*exp(-h/H))*v^2*Cd*Aref)^2)/(mass*g))-gmax)/gmax)','m/s','gpart'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'epsQ*sec(0.5*pi*(2*k*sqrt((rho0*exp(-h/H))/rn)*v^3-Qmax)/Qmax)','m/s','Qpart'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {'h/1000','km','Altitude'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {'k*sqrt(rho0*exp(-h/H)/rn)*v^3/100^2','W/cm^2','Heat rate'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'h/1000','km','Altitude'};
% in.figure(ind).plot(1).y{1} = {'0.5*rho0*exp(-h/H)*v^2','kg/(m*s^2)','Dynamic Pressure'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% rho = '(rho0*exp(-(r-re)/H))'; % Exponential Atmospheric Density [kg/m^3] 
% Cl = '(1.5658*alfa + -0.0000)'; 
% Cd = '(1.6537*alfa^2 + -0.0000*alfa + 0.0612)';
% D = ['(1/2*',rho,'*v^2*',Cd,'*Aref)']; % Drag Force [N] 
% L = ['(1/2*',rho,'*v^2*',Cl,'*Aref)']; % Lift Force [N] 
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {['sqrt(',L,'^2+',D,'^2)/mass/9.81'],'','G-Loading'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];

% %%%%%%%%%%%%%%
% %% Subplots %%
% %%%%%%%%%%%%%%
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'t','s','Time'};
% in.figure(ind).plot(1).y{1} = {'lamR/1000','km/s^2','\lambda_R'};
% in.figure(ind).plot(1).legend.location = 'NorthEast';
% in.figure(ind).plot(2).x{1} = {'t','s','Time'};
% in.figure(ind).plot(2).y{1} = {'lamTHETTA/1000000','km^2/s^2/rad','\lambda_\theta'};
% in.figure(ind).plot(2).legend.location = 'NorthWest';
% in.figure(ind).plot(3).x{1} = {'t','s','Time'};
% in.figure(ind).plot(3).y{1} = {'lamV/1000','km/s','\lambda_v'};
% in.figure(ind).plot(3).legend.location = 'SouthEast';
% in.figure(ind).plot(4).x{1} = {'t','s','Time'};
% in.figure(ind).plot(4).y{1} = {'lamGAM/1000000','km^2/s^2/rad','\lambda_\gamma'};
% in.figure(ind).plot(4).legend.location = 'NorthWest';
% in.figure(ind).plot(5).x{1} = {'t','s','Time'};
% in.figure(ind).plot(5).y{1} = {'lamPHI/1000000','km^2/s^2/rad','\lambda_\phi'};
% in.figure(ind).plot(5).legend.location = 'NorthEast';
% in.figure(ind).plot(6).x{1} = {'t','s','Time'};
% in.figure(ind).plot(6).y{1} = {'lamPSII/1000000','km^2/s^2/rad','\lambda_\psi'};
% in.figure(ind).plot(6).legend.location = 'NorthEast';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
% in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
% in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
% in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'alfa*180/pi','deg','Angle of Attack'};
% in.figure(ind).plot(1).legend.location = 'NorthEast';
% in.figure(ind).plot(2).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(2).y{1} = {'bank*180/pi','deg','Bank Angle'};
% in.figure(ind).plot(2).legend.location = 'NorthWest';
% % in.figure(ind).plot(3).x{1} = {'t',0,1,'s','Time'};
% % in.figure(ind).plot(3).y{1} = {'d2Hdu2',0,1/1000000,'km^2/s^2/rad^2','d^2H/du^2'};
% % in.figure(ind).plot(3).legend.location = 'NorthWest';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
% in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
% in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
% in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'(r-6378000)/1000','km','Altitude'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% in.figure(ind).plot(2).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(2).y{1} = {'bank*180/pi','deg','Bank Angle'};
% in.figure(ind).plot(2).legend.location = 'SouthWest';
% in.figure(ind).plot(3).x{1} = {'thetta*180/pi','deg','Longitude'};
% in.figure(ind).plot(3).y{1} = {'phi*180/pi','deg','Latitude'};
% in.figure(ind).plot(3).legend.location = 'NorthEast';
% in.figure(ind).plot(4).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(4).y{1} = {'alfa*180/pi','deg','Angle of Attack'};
% in.figure(ind).plot(4).legend.location = 'NorthEast';
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = 'testMovie';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
% in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
% in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
% in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];
% 
% ind = ind+1;
% in.figure(ind).plot(1).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(1).y{1} = {'(r-6378000)/1000','km','Altitude'};
% in.figure(ind).plot(1).legend.location = 'NorthWest';
% % in.figure(ind).plot(1).axis = [2 7 0 100];
% in.figure(ind).plot(2).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(2).y{1} = {'bank*180/pi','deg','Bank Angle'};
% in.figure(ind).plot(2).legend.location = 'SouthWest';
% % in.figure(ind).plot(2).axis = [2 7 -200 200];
% in.figure(ind).plot(3).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(3).y{1} = {'k*sqrt(rho0*exp(-(r-re)/H)/rn)*v^3/100^2','W/cm^2','Heat rate'};
% in.figure(ind).plot(3).legend.location = 'NorthEast';
% % in.figure(ind).plot(3).axis = [2 7 0 5000];
% in.figure(ind).plot(4).x{1} = {'v/1000','km/s','Velocity'};
% in.figure(ind).plot(4).y{1} = {'alfa*180/pi','deg','Angle of Attack'};
% in.figure(ind).plot(4).legend.location = 'NorthEast';
% % in.figure(ind).plot(4).axis = [2 7 0 25];
% in.figure(ind).movie.make = false;
% in.figure(ind).movie.name = '/home/mjgrant/outbox/movieHeatRateConstraint';
% in.figure(ind).movie.figPos = [1 31 1280 856*0.9];
% in.figure(ind).movie.plot(1).pos = [0.07 0.57 0.38 0.38];
% in.figure(ind).movie.plot(2).pos = [0.07 0.07 0.38 0.38];
% in.figure(ind).movie.plot(3).pos = [0.57 0.57 0.38 0.38];
% in.figure(ind).movie.plot(4).pos = [0.57 0.07 0.38 0.38];

return

