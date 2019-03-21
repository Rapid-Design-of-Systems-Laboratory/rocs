function [IG] = defaultGuess(in,varargin)
%
% This function calculates the initial guess for the easier optimal
% control problem that uses indirect method to evaluate the optimum answer.
%
% input : in [structures]
%         in is obtained from inputsRunTrajectoryProcess.m
% output : IG [structure]
% Developed by : Dr. M.J. Grant, Kshitij Mall and Thomas Antony

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%
numStates = length(in.oc.state(:,1)); % Total number of states
numControl = length(in.oc.control(:,1)); % Total number of controls

% Assume in.oc.guess.initial and in.oc.guess.terminal contains all the required values

% Initialize constraints on states (overwrites reverse integration states)
% in.oc.guess.initial.x = 0;
% in.oc.guess.initial.y = 0;
% in.oc.guess.initial.v = 1;
% in.oc.guess.initial.thetta = -pi/2;
%
% in.oc.guess.terminal.x =   in.oc.guess.initial.x + 0.1;
% in.oc.guess.terminal.y =   in.oc.guess.initial.y - 0.1;
% in.oc.guess.terminal.v = 2*in.oc.guess.initial.v + 0;
% in.oc.guess.terminal.thetta = 0;

numNuTotal = length(in.oc.constraint.initial(:,1)) + length(in.oc.constraint.terminal(:,1));

% in.oc.guess.timeIntegrate = 0.1;

% Set initial point
guessX = nan(2*numStates,1); % Initialize guess vector

% Set parameters of integration
if in.oc.guess.timeIntegrate > 0 % Forward integrate
  tSpan = [0 1]; % Normalized time
  % States
  for ctr = 1 : 1 : numStates
    guessX(ctr) = in.oc.guess.initial.([char(in.oc.state(ctr))]);
  end
else % Reverse integrate
  tSpan = [0 -1]; % Normalized time
  % States
  for ctr = 1 : 1 : numStates
    guessX(ctr) = in.oc.guess.terminal.([char(in.oc.state(ctr))]);
  end
end
options = in.odeOptions; % Options for ODE solver

% Costates
guessX(numStates+1:end,1) = in.oc.guess.costate*ones(numStates,1); 

% Arbitrarily set x0 and xf (only used in boundary condition function, but
% arguments passed to derivative function must be equivalent)
x0 = NaN(numStates,1);
xf = NaN(numStates,1);

p_ind = 1;
if ~in.convertParametersToStates
	p(p_ind,1) = abs(in.oc.guess.timeIntegrate);
	p_ind = p_ind + 1;
else
	guessX(2*numStates+1) = abs(in.oc.guess.timeIntegrate);
end
p(p_ind:numNuTotal+p_ind-1,1) = zeros(numNuTotal,1); % Initialize remaining elements of p

% Integrate trajectory
arcSequence = 0;
interiorPointConstraintSequence = [];
interiorPointNumLagrangeMultipliers = [];


% Save values from constants
if isfield(in,'const')
  names = fieldnames(in.const);
else
  names = {};
end

for ctr1 = 1 : 1 : length(names)
  const.(char(names{ctr1})) = in.const.(char(names{ctr1})){1};
end

% Save values from constraints
if isfield(in,'constraintVal')
  names = fieldnames(in.constraintVal);
else
  names = {};
end

constraintVal = struct; % intitialize empty structure
for ctr1 = 1 : 1 : length(names)
  constraintVal.(char(names{ctr1})) = in.constraintVal.(char(names{ctr1})){1};
end
    

% Note that the derivative function is derivFunc_mex
[tNormalizedTemp,Xtemp] = ode45(@derivFunc,tSpan,guessX, ...
        options,p,const,constraintVal,arcSequence,interiorPointConstraintSequence, ...
        interiorPointNumLagrangeMultipliers,x0,xf);
      
% Flip back to forward propagation (in case of reverse integration)
if in.oc.guess.timeIntegrate < 0
  tNormalized = flipud(tNormalizedTemp)-min(tNormalizedTemp);
  X = flipud(Xtemp);
else % No alteration for forward integration
  tNormalized = tNormalizedTemp;
  X = Xtemp;
end


ind = 0;
if ~in.convertParametersToStates
	ind = ind + 1; p(ind) = abs(in.oc.guess.timeIntegrate);
end
for pIndex = 1:length(in.oc.constraint.initial(:,1))
    ind = ind+1;
    % Fix this!!!
    p(ind) = X(1,numStates+pIndex);
end

for pIndex = 1:length(in.oc.constraint.terminal(:,1))
    ind = ind+1;
    p(ind) = X(end,numStates+pIndex); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Initial Guess Structure %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sol.x = tNormalized';
sol.y = X';

% Save values from constants
if isfield(in,'const')
	names = fieldnames(in.const);
else
	names = {};
end

for ctr1 = 1 : 1 : length(names)
	eval(['const.',char(names{ctr1}),' = in.const.',char(names{ctr1}),'{1};']);
end

% Save values from constraints
if isfield(in,'constraintVal')
	names = fieldnames(in.constraintVal);
else
	names = {};
end

constraint = struct; % intitialize empty structure
for ctr1 = 1 : 1 : length(names)
	eval(['constraint.',char(names{ctr1}),' = in.constraintVal.',char(names{ctr1}),'{1};']);
end

% Compute control history
sol.control = NaN(numControl,length(sol.x));
sol.d2Hdu2 = NaN(1,length(sol.x));
% sol.d2Hdu2 = NaN(numControl,length(sol.x));
controlOutput = cell(1,1+numControl);
numArcs = 1;

for timeCtr = 1 : 1 : length(sol.x)
	if ~in.convertParametersToStates
		[controlOutput{:}] = computeControlUnconstrained(sol.y(:,timeCtr),const,constraint,numArcs);
	else
		[controlOutput{:}] = computeControlUnconstrained([sol.y(:,timeCtr);zeros(in.maxNumArcs-1,1)],const,constraint,1);
	end
	for ctr1 = 1 : 1 : numControl
		sol.control(ctr1,timeCtr) = controlOutput{ctr1};
	end
	% sol.d2Hdu2(:,timeCtr) = controlOutput{numControl+1};
end

sol.parameters = p;

inTemp = in;

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

IG.sol = sol;
% IG.in = inTemp;

IG.numStates = numStates;

return