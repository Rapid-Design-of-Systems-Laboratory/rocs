function testHigherOrderSTMmatrix3
% Example on page 59 in Bryson and Ho

close all; clc;

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

strEOM{1,1} = 'x-lam'; % x+lam^2 not exact?
strEOM{2,1} = 'lam^(4/5)';    % with -lam?
states{1,1} = 'x';
states{2,1} = 'lam';
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
tSpan = [0 10];
x0nom = 1*ones(length(states),1);
orderSTT = 4;
dx0 = 0.2*ones(length(states),1);
plotSet = {'b','k','r','g','m','c'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Jacobian Matrices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create symbolic variables
numStates = length(strEOM);
f = sym(sym(strEOM));
x = sym(sym(states));

tempJac = f;
strJac = cell(orderSTT,1);
for ctr = 1 : 1 : orderSTT
  
  strJac{ctr} = sym([]);
  
  % Determine size of Jacobian values are being assigned to
  sizeJac = numStates*ones(1,ctr+1);
  
  % Construct empty cell array to store subscript indices
  sub = cell(1,length(sizeJac));
  
  % Compute Jacobian
  for ind = 1 : 1 : numStates^(ctr+1)
    
    % Construct cell array of subscript elements in Jacobian    
    [sub{:}] = ind2sub(sizeJac,ind);
    
    % Compute Jacobian
    strJac{ctr}(sub{:}) = diff(tempJac(sub{1:end-1}),x(sub{end}));
    
    % Construct initial values for STT. Identity for first order STT and zero 
    % for rest.
    if (ctr == 1) && (sum(diff([sub{:}])) == 0)
      val = 1;
    else
      val = 0;
    end
    X0{ctr}(ind,1) = val;
    
  end
  
  % Save Jacobian
  tempJac = strJac{ctr};
  
end

%%%%%%%%%%%%%%
%% Run ODEs %%
%%%%%%%%%%%%%%

% Propagate nominal solution
[tNom,xNom] = ode45(@eomExact,tSpan,x0nom,options,states,strEOM);

% Propagate perturbed solution
[tPert,xPert] = ode45(@eomExact,tNom,x0nom+dx0,options,states,strEOM);

% Construct STT initial guess
X0vec = [];
for ctr = 1 : 1 : orderSTT
  X0vec = [X0vec; X0{ctr}];
end

% Propagate STT solutions
[tSTT,xSTT] = ode45(@eomSTT,tSpan,X0vec,options,states,strJac,tNom,xNom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Approximate STT Solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cycle through all nominal trajectory times
for tCtr = 1 : 1 : length(tNom)
  
  % Interpolate STT values
  STTt = interp1(tSTT,xSTT,tNom(tCtr))';
  
  % Reshape STT for easy indexing
  startInd = 1;
  for ctr = 1 : 1 : orderSTT
    
    sizeSTT = numStates*ones(1,ctr+1);
    STT{ctr} = reshape(STTt(startInd:numStates^(ctr+1)+startInd-1,1),sizeSTT);
    startInd = numStates^(ctr+1)+startInd;
    
  end
  
  % Compute approximate solution
  for i = 1 : 1 : numStates
  for p = 1 : 1 : orderSTT
    
    % Determine size of Jacobian values are being assigned to
    sizeSTT = numStates*ones(1,p);
    
    % Construct empty cell array to store subscript indices
    sub = cell(1,length(sizeSTT));
    
    % Compute change in solution for each term in STT
    sumVal = 0;
    for ind = 1 : 1 : numStates^p
      if p == 1
        [sub{:}] = ind2sub([sizeSTT 1],ind); % Find subscripts for STT
      else
        [sub{:}] = ind2sub(sizeSTT,ind); % Find subscripts for STT
      end
      subSTT = {i sub{:}}; %#ok<CCAT> % Add state to subscripts
      sumVal = sumVal + STT{p}(subSTT{:})*prod(dx0([sub{:}]));
    end
    sumVal = sumVal/factorial(p);
    dx{p}(tCtr,i) = sumVal;
    
  end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute STT-estimated Solution and Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ctr = 1 : 1 : numStates
  xPS_sum = xNom(:,ctr);
  for p = 1 : 1 : orderSTT
    xPS_sum = xPS_sum + dx{p}(:,ctr);
    xPS{p}(:,ctr) = xPS_sum;
    errPS{p}(:,ctr) = xPS{p}(:,ctr) - xPert(:,ctr);
  end
end

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

% Plot solution
figure(1);

% Plot nominal solution
for ctr = 1 : 1 : numStates
  subplot(numStates,1,ctr);
  hold on; grid on;
  plot(tNom,xNom(:,ctr),plotSet{1});
  xlabel('t');
  ylabel(states{ctr,1});
  title([states{ctr,1},' vs. t']);
end

% Plot perturbed solution
for ctr = 1 : 1 : numStates
  subplot(numStates,1,ctr);
  plot(tPert,xPert(:,ctr),plotSet{2});
  xlabel('t');
  ylabel(states{ctr,1});
  title([states{ctr,1},' vs. t']);
end

% Plot STT-estimated solution
for ctr = 1 : 1 : numStates
  subplot(numStates,1,ctr);
  for p = 1 : 1 : orderSTT
    plot(tNom,xPS{p}(:,ctr),plotSet{p+2});
  end
end

% Plot error
figure(2);
for ctr = 1 : 1 : numStates
  subplot(numStates,1,ctr);
  hold on; grid on;
  xlabel('t');
  ylabel(['Error in ',states{ctr,1}]);
  title(['Error in ',states{ctr,1},' vs. t']);
  for p = 1 : 1 : orderSTT
    plot(tNom,errPS{p}(:,ctr),plotSet{p+2});
  end
end

% Plot discrepancy between approximate solutions
figure(3);
for ctr = 1 : 1 : numStates
  subplot(numStates,1,ctr);
  hold on; grid on;
  xlabel('t');
  ylabel(['Discrepancy in ',states{ctr,1}]);
  title(['Discrepancy in ',states{ctr,1},' vs. t']);
  for p = 1 : 1 : orderSTT-1
    plot(tNom,xPS{p+1}(:,ctr)-xPS{p}(:,ctr),plotSet{p+2});
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xDot] = eomExact(t,xVec,states,strEOM) %#ok<INUSL>

% Inputs
numStates = length(states);

% Assign states
for ctr = 1 : 1 : numStates
  eval([states{ctr,1},'=xVec(ctr);']);
end

% Compute derivatives
for ctr = 1 : 1 : numStates
  eval(['xDot(ctr,1)=',strEOM{ctr,1},';']);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xDot] = eomSTT(t,xVec,states,strJac,tNom,xNom) %#ok<INUSD,INUSL>

% Inputs
numStates = length(states);
orderSTT = length(strJac);
% Interpolate nominal solution at time t
for ctr = 1 : 1 : numStates
  eval([states{ctr,1},'=interp1(tNom,xNom(:,ctr),t);']);
end

% Evaluate Jacobians which are only a function of the states and reshape input 
% from column vector to matrices for easy indexing in computations
startInd = 1;
for ctr = 1 : 1 : orderSTT
  
  % Determine size of Jacobian values are being assigned to
  sizeJac = numStates*ones(1,ctr+1);
  
  % Construct empty cell array to store subscript indices
  sub = cell(1,length(sizeJac));
  
  % Compute Jacobian
  for ind = 1 : 1 : numStates^(ctr+1)
    
    % Construct cell array of subscript elements in Jacobian    
    [sub{:}] = ind2sub(sizeJac,ind);
    
    % Compute Jacobian
    valJac{ctr}(sub{:}) = eval(strJac{ctr}(sub{:}));
    
  end
  
  % Construct STT into matrix form for easy indexing
  dimArg = num2cell(sizeJac);
  STT{ctr} = reshape(xVec(startInd:numStates^(ctr+1)+startInd-1,1),dimArg{:});
    
  % Increment start index
  startInd = numStates^(ctr+1)+startInd;
  
end

% Compute STT rates
for ctr = 1 : 1 : orderSTT
  if ctr == 1
    for i = 1 : 1 : numStates
    for a = 1 : 1 : numStates
      % Perform summation for term
      sumVal = 0;
      for alpha = 1 : 1 : numStates
        sumVal = sumVal + valJac{ctr}(i,alpha)*STT{ctr}(alpha,a);
      end
      STTdot{ctr}(i,a) = sumVal;
    end
    end
  elseif ctr == 2
    for i = 1 : 1 : numStates
    for a = 1 : 1 : numStates
    for b = 1 : 1 : numStates
      % Perform summation for term
      sumVal = 0;
      for alpha = 1 : 1 : numStates
        sumVal = sumVal + valJac{ctr-1}(i,alpha)*STT{ctr}(alpha,a,b);
      for beta = 1 : 1 : numStates
        sumVal = sumVal + ...
          valJac{ctr}(i,alpha,beta)*STT{ctr-1}(alpha,a)*STT{ctr-1}(beta,b);
      end
      end
      STTdot{ctr}(i,a,b) = sumVal;
    end
    end
    end
  elseif ctr == 3
    for i = 1 : 1 : numStates
    for a = 1 : 1 : numStates
    for b = 1 : 1 : numStates
    for c = 1 : 1 : numStates
      % Perform summation for term
      sumVal = 0;
      for alpha = 1 : 1 : numStates
        sumVal = sumVal + valJac{ctr-2}(i,alpha)*STT{ctr}(alpha,a,b,c);
      for beta = 1 : 1 : numStates
%         % Old way
%         sumVal = sumVal + valJac{ctr-1}(i,alpha,beta)* ...
%           (STT{ctr-2}(alpha,a)*STT{ctr-1}(beta,b,c) + ...
%           STT{ctr-1}(alpha,a,b)*STT{ctr-2}(beta,c) + ...
%           STT{ctr-1}(alpha,a,c)*STT{ctr-2}(beta,b));
        sumVal = sumVal + valJac{ctr-1}(i,alpha,beta)* ...
          3/2*(STT{ctr-2}(alpha,a)*STT{ctr-1}(beta,b,c) + ...
          STT{ctr-1}(alpha,a,b)*STT{ctr-2}(beta,c));
      for gamma = 1 : 1 : numStates
        sumVal = sumVal + ...
          valJac{ctr}(i,alpha,beta,gamma)*STT{ctr-2}(alpha,a)* ...
          STT{ctr-2}(beta,b)*STT{ctr-2}(gamma,c);
      end
      end
      end
      STTdot{ctr}(i,a,b,c) = sumVal;
    end
    end
    end
    end
  elseif ctr == 4
    for i = 1 : 1 : numStates
    for a = 1 : 1 : numStates
    for b = 1 : 1 : numStates
    for c = 1 : 1 : numStates
    for d = 1 : 1 : numStates
      % Perform summation for term
      sumVal = 0;
      for alpha = 1 : 1 : numStates
        sumVal = sumVal + valJac{ctr-3}(i,alpha)*STT{ctr}(alpha,a,b,c,d);
      for beta = 1 : 1 : numStates
%         % Old way
%         sumVal = sumVal + valJac{ctr-2}(i,alpha,beta)* ...
%           (STT{ctr-1}(alpha,a,b,c)*STT{ctr-3}(beta,d) + ...
%           STT{ctr-1}(alpha,a,b,d)*STT{ctr-3}(beta,c) + ...
%           STT{ctr-1}(alpha,a,c,d)*STT{ctr-3}(beta,b) + ...
%           STT{ctr-2}(alpha,a,b)*STT{ctr-2}(beta,c,d) + ...
%           STT{ctr-2}(alpha,a,c)*STT{ctr-2}(beta,b,d) + ...
%           STT{ctr-2}(alpha,a,d)*STT{ctr-2}(beta,b,c) + ...
%           STT{ctr-3}(alpha,a)*STT{ctr-1}(beta,b,c,d));
          sumVal = sumVal + valJac{ctr-2}(i,alpha,beta)* ...
            (2*(STT{ctr-3}(alpha,a)*STT{ctr-1}(beta,b,c,d) + ...
            STT{ctr-3}(beta,a)*STT{ctr-1}(alpha,b,c,d)) + ...
            3*(STT{ctr-2}(alpha,a,b)*STT{ctr-2}(beta,c,d)));
      for gamma = 1 : 1 : numStates
%         % Old way
%         sumVal = sumVal + valJac{ctr-1}(i,alpha,beta,gamma)* ...
%           (STT{ctr-2}(alpha,a,b)*STT{ctr-3}(beta,c)*STT{ctr-3}(gamma,d) + ...
%           STT{ctr-2}(alpha,a,c)*STT{ctr-3}(beta,b)*STT{ctr-3}(gamma,d) + ...
%           STT{ctr-2}(alpha,a,d)*STT{ctr-3}(beta,b)*STT{ctr-3}(gamma,c) + ...
%           STT{ctr-3}(alpha,a)*STT{ctr-2}(beta,b,c)*STT{ctr-3}(gamma,d) + ...
%           STT{ctr-3}(alpha,a)*STT{ctr-2}(beta,b,d)*STT{ctr-3}(gamma,c) + ...
%           STT{ctr-3}(alpha,a)*STT{ctr-3}(beta,b)*STT{ctr-2}(gamma,c,d));
        sumVal = sumVal + valJac{ctr-1}(i,alpha,beta,gamma)* ...
          2*(STT{ctr-3}(alpha,a)*STT{ctr-3}(beta,b)*STT{ctr-2}(gamma,c,d) + ...
          STT{ctr-3}(alpha,a)*STT{ctr-2}(beta,b,c)*STT{ctr-3}(gamma,d) + ...
          STT{ctr-2}(alpha,a,b)*STT{ctr-3}(beta,c)*STT{ctr-3}(gamma,d));
      for delta = 1 : 1 : numStates
        sumVal = sumVal + valJac{ctr}(i,alpha,beta,gamma,delta)* ...
          STT{ctr-3}(alpha,a)*STT{ctr-3}(beta,b)* ...
          STT{ctr-3}(gamma,c)*STT{ctr-3}(delta,d);
      end
      end
      end
      end
      STTdot{ctr}(i,a,b,c,d) = sumVal;
    end
    end
    end
    end
    end
  end
end

% Reshape STT rates to column vector for output
startInd = 1;
for ctr = 1 : 1 : orderSTT
  numRows = numStates^(ctr+1);
  xDot(startInd:numStates^(ctr+1)+startInd-1,1) = ...
    reshape(STTdot{ctr},numRows,1);
  startInd = startInd + numStates^(ctr+1);
end

return

