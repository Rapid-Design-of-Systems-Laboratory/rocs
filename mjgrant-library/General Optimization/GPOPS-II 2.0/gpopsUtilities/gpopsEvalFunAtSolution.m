function [contoutput, endpoutput] = gpopsEvalFunAtSolution(result)

% ---------------------------------------------------------------- %
% GPOPSEVALUSERFUN                                                 %
% This function evaluates the continuous and endpoint functions on % 
% the solution.                                                    %
% ---------------------------------------------------------------- %

% ---------------------------%
% Extract Setup and Solution %
% ---------------------------%
setup = result.setup;
solution = result.solution;

% -------------------- %
% Get Number of Phases %
% -------------------- %
numphase = length(solution.phase);

% ------------------------ %
% Get Number of Parameters %
% ------------------------ %
if isfield(solution,'parameter');
    numparameter = length(solution.parameter);
else
    numparameter = 0;
end

% ------------------------------------ %
% Pre-allocate CONTPHASE and ENDPPHASE % 
% ------------------------------------ %
contphase(numphase).state = [];
endpphase(numphase).initialstate = [];
endpphase(numphase).finalstate = [];
if isfield(solution.phase, 'control');
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
if isfield(solution.phase, 'integral');
    integralswitch = true;
    endpphase(numphase).integral = [];
else
    integralswitch = false;
end

% ------------------------------- %
% Get Input Values for Each Phase %
% ------------------------------- %
for phasecount = 1:numphase;
  % --------------------------- %
  % Get Solution for Each Phase %
  % --------------------------- %
  solphase = solution.phase(phasecount);
  
  % ---------------------------- %
  % Get Input for State Solution %
  % ---------------------------- %
  contphase(phasecount).state = solphase.state;
  endpphase(phasecount).initialstate = solphase.state(1,:);
  endpphase(phasecount).finalstate = solphase.state(end,:);
  
  % ------------------------------ %
  % Get Input for Control Solution %
  % ------------------------------ %
  if controlswitch;
    contphase(phasecount).control = solphase.control;
  end
  
  % --------------------------- %
  % Get input for Time Solution %
  % --------------------------- %
  contphase(phasecount).time = solphase.time;
  endpphase(phasecount).initialtime = solphase.time(1);
  endpphase(phasecount).finaltime = solphase.time(end);
  
  % ------------------------------------ %
  % Get CONTPHASE for PARAMETER Solution %
  % ------------------------------------ %
  if numparameter ~= 0;
    contphase(phasecount).parameter = ones(length(solphase.time),1)*solution.parameter;
  end
  
  % ----------------------------------- %
  % Get ENDPPHASE for integral solution %
  % ----------------------------------- %
  if integralswitch;
    endpphase(phasecount).integral = solphase.integral;
  end
end

% --------------------------------------------------------- %
% Adding Solution for All Phases to CONTINPUT and ENDPINPUT %
% --------------------------------------------------------- %
continput.phase = contphase;
endpinput.phase = endpphase;

% ------------------------------------ %
% Get ENDPINPUT for parameter solution %
% ------------------------------------ %
if numparameter ~= 0;
    endpinput.parameter = solution.parameter;
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