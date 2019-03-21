function [in] = spaceShuttle_hc() 
% 
% This function creates input for the optimization problem,  
% which is later used by other functions  
% 
% input : void 
% output : in [structure] 
% Developed by : Dr. M.J. Grant and Kshitij Mall 
% Last modified: Jan 20, 2017 

%%%%%%%%%%%%%%%%%%%%%%% 
%% Execution Control %% 
%%%%%%%%%%%%%%%%%%%%%%% 

in.bvpOptions = bvpset('AbsTol',1e-4,'RelTol',1e-4,'Nmax',100000000,'Stats','on'); 
in.oc.writeEquations = false;

%%%%%%%%%%%%%
%% Scaling %%
%%%%%%%%%%%%%

% Doubles must be equal to one (not dynamically updated during continuation)
in.autoScale = true;
in.scale = {'m','x1'; ...
		    'rad',1; ...
			's','x1/x3'; ...
			'nd',1}; % nd = nondimensional

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Independent Variable %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g., time

in.oc.independentVariable = {'t','s'}; % time

%%%%%%%%%%%%
%% States %%
%%%%%%%%%%%%

in.oc.state = {'h','m'; ... % radial position magnitude
			   'phi','rad'; ... % latitude, positive northward
			   'v','m/s'; ... % relative velocity
			   'gam','rad'; ... % relative flight-path angle
			   'psii','rad'}; % azimuth

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of Motion %%
%%%%%%%%%%%%%%%%%%%%%%%%%
hbar = '(h/hs - 1)';
H1 = '(b^2*h^2/v^2)';
H2 = ['(b*h/v - ',H1,')'];
H3 = ['(1 - b*h/v - ',H2,')'];
H4 = ['(v/(b*h) - 2 + b*h/v - ',H3,')'];
H5 = ['(v^2/(b*h)^2 - 3*v/(b*h) + 3 - b*h/v -',H4,')'];
B1 = ['(0.110717 + 0.834519*',hbar,' + 1.213679*',hbar,'^2 - 1.060833*',hbar,'^3)'];
B2 = ['(-0.672677 + 2.73417*',hbar,' - 0.864369*',hbar,'^2 - 12.1*',hbar,'^3)'];
B3 = ['(0.812241 + 2.337815*',hbar,' + 10.31628*',hbar,'^2 + 22.97486*',hbar,'^3)'];
B4 = ['(-3.151267 - 13.62131*',hbar,' - 40.4855*',hbar,'^2 - 57.83333*',hbar,'^3)'];
B5 = ['(2.368095 + 19.0734*',hbar,' + 69.86905*',hbar,'^2 + 127.777778*',hbar,'^3)'];
CLH = ['(',B1,'*',H1,'+',B2,'*',H2,'+',B3,'*',H3,'+',B4,'*',H4,'+',B5,'*',H5,'+ delCLH)']

CL = ['(0.5*(((1-w)*CLUB+w*',CLH,'-CLLB)*sin(CLw)+ (1-w)*CLUB+w*',CLH,'+CLLB))'];
Cd = ['(',CL,'^1.86 + 0.04)'];

qA = '(SA*exp(-h/H)*v^2)'; % Exponential atmospheric density
r = '(re+h)'; % Radial distance
D = ['(',qA,'*',Cd,')']; % Drag Force
L = ['(',qA,'*',CL,')']; % Lift Force

% 3D
Ft = D; % force along velocity vector
Fn = L; % force perpendicular to velocity vector
bank = '(pi*(1+sin(bankt))/4)'; 
in.oc.stateRate = {'v*sin(gam)'; ...
				   ['v*cos(gam)*sin(psii)/',r]; ...
				   ['-',Ft,'- mu*sin(gam)/',r,'^2']; ... 
  				   [Fn,'*cos(',bank,')/v  - mu/(v*',r,'^2)*cos(gam) + v/',r,'*cos(gam)']; ...
  				   [Fn,'*sin(',bank,')/(cos(gam)*v) - v/',r,'*cos(gam)*cos(psii)*tan(phi)']};

%%%%%%%%%%%%%
%% Control %%
%%%%%%%%%%%%%

in.oc.control = {'bankt','rad';... % bank angle control
				 'CLw','nd'}; % angle of attack control
in.oc.assumptions.control = ''; % assumptions for Mathematica when solving for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path cost
in.oc.cost.path = {'0','rad'}; 

% Terminal cost 
in.oc.cost.terminal = {'-phi','rad'}; 

% Initial cost 
in.oc.cost.initial = {'0','rad'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Endpoint Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial constraint 
in.oc.constraint.initial = {'h-x0(1)','m'; ...
							'phi-x0(2)','rad'; ...
							'v-x0(3)','m/s'; ...
							'gam-x0(4)','rad';...
							'psii-x0(5)','rad'};

% Final constraint 
in.oc.constraint.terminal = {'h-xf(1)','m'; ...
							 'v-xf(3)','m/s'; ...
							 'gam-xf(4)','rad'};                       

%%%%%%%%%%%%%%%
%% Constants %%
%%%%%%%%%%%%%%%

in.const.mu = {9.80665*6371200^2,'m^3/s^2'}; % Gravitational parameter
in.const.H = {1000/0.145,'m'}; % Scale height for atmosphere of Earth
in.const.re = {6371200,'m'}; % Radius of planet
in.const.SA = {3.08/1000,'1/m'}; % Reference area of vehicle
in.const.CLLB = {0.1,'nd'}; % Reference area of vehicle
in.const.CLUB = {0.3,'nd'}; % Reference area of vehicle
in.const.delCLH = {0.12,'nd'}; % Reference area of vehicle
in.const.hs = {50000,'m'}; % Reference area of vehicle
in.const.b = {0.095,'1/s'}; % Reference area of vehicle
in.const.w = {0,'nd'}; % Reference area of vehicle

% Ham = ['lamH*(v*sin(gam)) + lamPHI*(v*cos(gam)*sin(psii)/',r,') + lamV*(-',Ft,'- mu*sin(gam)/',r,'^2) + lamGAM*(',Fn,'*cos(',bank,')/v  - mu/(v*',r,'^2)*cos(gam) + v/',r,'*cos(gam)) + lamPSII*(',Fn,'*sin(',bank,')/(cos(gam)*v) - v/',r,'*cos(gam)*cos(psii)*tan(phi))'];

%%%%%%%%%%%%%%%%%%%
%% Initial Guess %%
%%%%%%%%%%%%%%%%%%%
in.oc.guess.mode = 'auto';

in.oc.guess.timeIntegrate = -50; % 0.1 leads to local min
% Terminal
in.oc.guess.terminal.h = 30000;
in.oc.guess.terminal.phi = 44.42*pi/180;
in.oc.guess.terminal.v = 1116;
in.oc.guess.terminal.gam = -2.7*pi/180;
in.oc.guess.terminal.psii = 81.19*pi/180;
in.oc.costates = -0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Execution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in.cont.method = 1; % 1 = manually changing parameters, 2 = targeting new solution using STTs  

ind = 0;

%%%%%%%%%%%%%%%%%%%%%%
%% Continuation Set %%
%%%%%%%%%%%%%%%%%%%%%%
% 
ind = ind+1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 100; % 41; Number of steps in the continuation set
in.CONT{ind}.constraint.initial.v = 7850;
in.CONT{ind}.constraint.initial.h = 95000;
in.CONT{ind}.constraint.initial.gam = -1.25*pi/180;
in.CONT{ind}.constraint.initial.psii = 0*pi/180;
in.CONT{ind}.constraint.initial.phi = 0*pi/180;

ind = ind + 1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 2; % Number of steps in the continuation set
in.CONT{ind}.const.w = linspace(0,1,in.CONT{ind}.numCases); 

ind = ind + 1;
in.cont.method(ind) = 1;
in.CONT{ind}.numCases = 150; % Number of steps in the continuation set
in.CONT{ind}.const.delCLH = linspace(0,-0.12-0.04566259,in.CONT{ind}.numCases); % -0.05
return
