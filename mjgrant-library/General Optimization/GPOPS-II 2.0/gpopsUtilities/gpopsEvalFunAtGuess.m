function [contoutput, endpoutput] = gpopsEvalFunAtGuess(setup)

% ------------------------------------------------------------- %
% GPOPSEVALFUNATGUESS                                           %
% This function evaluates the continuous and endpoint functions %
% on the guess.                                                 %
% ------------------------------------------------------------- %


% ---------------- %
%    Get guess     %
% ---------------- %
guess = setup.guess;

% ------------------------ %
% Get the Number of Phases %
% ------------------------ %
numphase = length(guess.phase);

% ------------------------ %
% Get Number of Parameters %
% ------------------------ %
if isfield(guess,'parameter');
    numparameter = length(guess.parameter);
else
    numparameter = 0;
end

% ------------------------------------ %
% Pre-allocate CONTPHASE and ENDPPHASE %
% ------------------------------------ %
contphase(numphase).state = [];
endpphase(numphase).initialstate = [];
endpphase(numphase).finalstate = [];
if isfield(guess.phase, 'control');
    controlswitch = true;
    contphase(numphase).control = [];
else
    controlswitch = false;
end
contphase(numphase).time = [];
endpphase(numphase).initialtime = [];
endpphase(numphase).finaltime = [];
if numparameter ~= 0;
    contphase(numphase).parameter = [];
end
if isfield(guess.phase, 'integral');
    integralswitch = true;
    endpphase(numphase).integral = [];
else
    integralswitch = false;
end

% ------------------------------------- %
%   Get the Input Values for Each Phase %
% ------------------------------------- %
for phasecount = 1:numphase;
  % ------------------------ %
  % Get Guess for Each Phase %
  % ------------------------ %
  guessphase = guess.phase(phasecount);
  
  % ------------------------- %
  % Get Input for State Guess %
  % ------------------------- %
  contphase(phasecount).state = guessphase.state;
  endpphase(phasecount).initialstate = guessphase.state(1,:);
  endpphase(phasecount).finalstate = guessphase.state(end,:);
  
  % --------------------------- %
  % Get Input for Control Guess %
  % --------------------------- %
  if controlswitch;
    contphase(phasecount).control = guessphase.control;
  end
  
  % ------------------------ %
  % Get Input for Time Guess %
  % ------------------------ %
  contphase(phasecount).time = guessphase.time;
  endpphase(phasecount).initialtime = guessphase.time(1);
  endpphase(phasecount).finaltime = guessphase.time(end);
  
  % --------------------------------- %
  % Get CONTPHASE for Parameter Guess %
  % --------------------------------- %
  if numparameter ~= 0;
    contphase(phasecount).parameter = ones(length(guessphase.time),1)*guess.parameter;
  end
  
  % -------------------------------- %
  % Get ENDPPHASE for Integral Guess %
  % -------------------------------- %
  if integralswitch;
    endpphase(phasecount).integral = guessphase.integral;
  end
end

% ------------------------------------------------------ %
% Adding Guess for All Phases to CONTINPUT and ENDPINPUT %
% ------------------------------------------------------ %
continput.phase = contphase;
endpinput.phase = endpphase;

% --------------------------------- %
% Get ENDPINPUT for Parameter Guess %
% --------------------------------- %
if numparameter ~= 0;
    endpinput.parameter = guess.parameter;
end

% -------------------------------------- %
% Add AUXDATA to CONTINPUT and ENDPINPUT %
% -------------------------------------- %
if isfield(setup,'auxdata');
    continput.auxdata = setup.auxdata;
    endpinput.auxdata = setup.auxdata;
end

% ---------------------------------------------------- %
% Evaluate Optimal Control Problem Continuous Function %
% ---------------------------------------------------- %
contoutput = feval(setup.functions.continuous, continput);

% -------------------------------------------------- %
% Evaluate Optimal Control Problem Endpoint Function %
% -------------------------------------------------- %
endpoutput = feval(setup.functions.endpoint, endpinput);
