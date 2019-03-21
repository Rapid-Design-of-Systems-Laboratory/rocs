function [IG] = getInitialGuess(in, oc)
%
% This function calculates the initial guess for the easier optimal
% control problem that uses indirect method to evaluate the optimum answer.
%
% input : in [structures]
%         in is obtained from inputsRunTrajectoryProcess.m
% output : IG [structure]
% Developed by : Dr. M.J. Grant and Kshitij Mall

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

numStates = length(in.oc.state(:,1)); % Total number of states
numControl = length(in.oc.control(:,1)); % Total number of controls

% For initial trajectory maximize terminal altitude
in.traj.h = 0; % Initial Radius 
in.traj.v = 0; % Initial Velocity 
in.traj.mass = 3;
in.ID.timeIntegrate = 10; 

% % For initial trajectory maximize terminal altitude
% in.traj.x1 = 0; % Initial Radius 
% in.traj.x2 = 0; % Initial Velocity 
% in.ID.timeIntegrate = -0.5; 

numNuTotal = 4; % For free radial magnitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Initial Guess Through Integration %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters of integration
if in.ID.timeIntegrate > 0 % Forward integrate
  tSpan = [0 1]; % Normalized time
else % Reverse integrate
  tSpan = [0 -1]; % Normalized time
end
options = in.odeOptions; % Options for ODE solver
p(1,1) = abs(in.ID.timeIntegrate);
p(2:numNuTotal+1,1) = zeros(numNuTotal,1); % Initialize remaining elements of p

% Set final point
guessX = nan(2*numStates,1); % Initialize guess vector

% States
for ctr = 1 : 1 : numStates
  
  guessX(ctr) = in.traj.([char(in.oc.state(ctr,1))]);
  
end

% Costates

% guessX(numStates+1:end,1) = zeros(numStates,1); 
guessX(numStates+1,1) = 0.1; 
guessX(numStates+2,1) = -0.1; 
guessX(numStates+3,1) = -0.1;
% Arbitrarily set x0 and xf (only used in boundary condition function, but
% arguments passed to derivative function must be equivalent)
x0 = nan(numStates,1);
xf = nan(numStates,1);

% RECENT ADDITION
arcSequence = [0];
interiorPointConstraintSequence = nan(1,1);
interiorPointNumLagrangeMultipliers = nan(1,1);

% Save values from constants
names = fieldnames(in.const);
% constArray = [];
for ctr1 = 1 : 1 : length(names)
	eval(['const.',char(names{ctr1}),' = in.const.',char(names{ctr1}),'{1};']);
end

% Save values from constraints
names = fieldnames(in.constraintVal);
constraint = struct; % intitialize empty structure
for ctr1 = 1 : 1 : length(names)
	eval(['constraint.',char(names{ctr1}),' = in.constraintVal.',char(names{ctr1}),'{1};']);
end
% Integrate trajectory
if in.convertParametersToStates
    guessX = [guessX; p(1,1)];
    p = p(2:end,1);
end
[tNormalizedTemp,Xtemp] = ode45(@derivFuncRegion,tSpan,guessX, ...
        options,1,p,const,constraint, ...
        arcSequence,interiorPointConstraintSequence,interiorPointNumLagrangeMultipliers,x0,xf);

% Flip back to forward propagation (in case of reverse integration)
if in.ID.timeIntegrate < 0
  tNormalized = flipud(tNormalizedTemp)-min(tNormalizedTemp);
  X = flipud(Xtemp);
else % No alteration for forward integration
  tNormalized = tNormalizedTemp;
  X = Xtemp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Initial Guess Structure %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sol.x = tNormalized';
sol.y = X';
sol.control = NaN(numControl,length(sol.x));
sol.d2Hdu2 = NaN(1,length(sol.x));
sol.parameters = p;

inTemp = in;

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

IG.sol = sol;
IG.in = inTemp;

if ~in.convertParametersToStates
	IG.numStates = numStates;
else
	% Total states = numStates + numNuTotal + 1 (for tf)
	IG.numStates = numStates + 1;
end

return
