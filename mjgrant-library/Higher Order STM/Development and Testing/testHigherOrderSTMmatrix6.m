function testHigherOrderSTMmatrix6
% Change in final time testing. All other previous files tested STT equations
% for fixed final time.

% Example on page 59 in Bryson and Ho.

close all; clc;

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

strEOM{1,1} = 'x^(4/5)*lam^(4/5)'; % x+lam^2 not exact?
strEOM{2,1} = 'lam^(4/5)';    % with -lam?
states{1,1} = 'x';
states{2,1} = 'lam';
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
tSpan = [0 1];
x0nom = 1*ones(length(states),1);
orderSTT = 2;
dx0 = 0.2*ones(length(states),1);
plotSet = {'b','k','r','g','m','c'};

orderXdot = 1; % Order of Taylor Series approximation for nth derivative of x
xDotCount = 2; % Number of time derivatives of x to estimate perturbed solution with different final time
dtf = 0.1; % For propagation of perturbed solution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Jacobian Matrices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is for STT computations. Save in own structure to keep code separate.

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

% Compute Jacobian again to order requested for dxDot calculations. Somewhat
% redundant with the above. Keep in own structure for easy coding for now.

% Create symbolic variables
f = sym(sym(strEOM));
x = sym(sym(states));

tempJac = f;
strJac2 = cell(orderXdot,1);
for ctr = 1 : 1 : orderXdot
  
  strJac2{ctr} = sym([]);
  
  % Determine size of Jacobian values are being assigned to
  sizeJac = numStates*ones(1,ctr+1);
  
  % Construct empty cell array to store subscript indices
  sub = cell(1,length(sizeJac));
  
  % Compute Jacobian
  for ind = 1 : 1 : numStates^(ctr+1)
    
    % Construct cell array of subscript elements in Jacobian    
    [sub{:}] = ind2sub(sizeJac,ind);
    
    % Compute Jacobian
    strJac2{ctr}(sub{:}) = diff(tempJac(sub{1:end-1}),x(sub{end}));
    
  end
  
  % Save Jacobian
  tempJac = strJac2{ctr};
  
end

% For each Jacobian associated with each order p, compute derivatives with
% respect to x for nth time derivative of dx
strJacXdot = cell(orderXdot,xDotCount);
for ctr = 1 : 1 : orderXdot
  
  strJacXdot{ctr,1} = strJac2{ctr};
  tempJac = strJacXdot{ctr,1};
  
  for ctr2 = 2 : 1 : xDotCount % No derivative on f for xDotCount = 1
  
    strJacXdot{ctr,ctr2} = sym([]);

    % Determine size of Jacobian values are being assigned to. (ctr+1)
    % corresponds to the size of the STT Jacobian. (ctr2-1) corresponds to the
    % number of time derivatives applied to the STT Jacobian.
    sizeJac = numStates*ones(1,(ctr+1)+(ctr2-1));

    % Construct empty cell array to store subscript indices
    sub = cell(1,length(sizeJac));

    % Compute Jacobian
    for ind = 1 : 1 : numStates^((ctr+1)+(ctr2-1))

      % Construct cell array of subscript elements in Jacobian    
      [sub{:}] = ind2sub(sizeJac,ind);

      % Compute Jacobian
      strJacXdot{ctr,ctr2}(sub{:}) = diff(tempJac(sub{1:end-1}),x(sub{end}));

    end
    
    tempJac = strJacXdot{ctr,ctr2};
  
  end
  
end

%%%%%%%%%%%%%%
%% Run ODEs %%
%%%%%%%%%%%%%%

% Propagate nominal solution
[tNom,xNom] = ode45(@eomExact,tSpan,x0nom,options,states,strEOM);

% Propagate perturbed solution
tSpanPert = [tSpan(1) tSpan(2)+dtf];
[tPert,xPert] = ode45(@eomExact,tSpanPert,x0nom+dx0,options,states,strEOM);

% Construct STT initial guess
X0vec = [];
for ctr = 1 : 1 : orderSTT
  X0vec = [X0vec; X0{ctr}];
end

% Obtain Taylor Series information for STTs
[taylorSet] = eomTermsSTT(orderSTT);

% Propagate nominal solution to perturbed trajectory time

% Estimate perturbed solution beyond nominal solution final time
[tfPertEst,xfPertEst] = dxfEstimate(strJacXdot,xDotCount,orderXdot);

% Propagate STT solutions
[tSTT,xSTT] = ode45(@eomSTT,tSpan,X0vec,options,states,strJac,taylorSet, ...
  tNom,xNom);

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
    % Interpolate error since trajectories vary in length
    xPertInterp = interp1(tPert,xPert(:,ctr),tNom);
    errPS{p}(:,ctr) = xPS{p}(:,ctr) - xPertInterp;
%     errPS{p}(:,ctr) = xPS{p}(:,ctr) - xPert(:,ctr);
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

function [xDot] = eomSTT(t,xVec,states,strJac,taylorSet,tNom,xNom) %#ok<INUSD,INUSL>

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

% Compute STT rates for each order in Taylor Series approximation
for p = 1 : 1 : orderSTT
  
  % Determine size and subscripts of STT
  sizeSTT = numStates*ones(1,p+1);
  subSTT = cell(1,length(sizeSTT));
  
  % Loop through each term in STT
  for indSTT = 1 : 1 : numStates^(p+1)
    
    % Determine subscript associated with STT
    [subSTT{:}] = ind2sub(sizeSTT,indSTT);
    
    % Loop through each Jacobian term. Perform summation for term associated 
    % with each Jacobian.
    sumVal = 0;
    for pJac = 1 : 1 : p

      % Determine number of terms asocciated with Jacobian
      numTerms = length(taylorSet{p,pJac}.term);
        
      % Determine size and subscripts associated with summation
      sizeSum = numStates*ones(1,pJac);
      subSum = cell(1,length(sizeSum));

      % Perform summation
      for indSum = 1 : 1 : numStates^pJac

        % Determine subscript associated with summation
        [subSum{:}] = ind2sub(sizeSum,indSum);

        % Add terms within summation
        sumTerms = 0;
        for ctrTerm = 1 : 1 : numTerms

          % Get order for each product in term
          term = taylorSet{p,pJac}.term(ctrTerm);

          % Determine number of products in each term
          numProd = length(term{:});

          % Perform multiplication of products
          prod = 1;
          termInd0 = 2; % Shift by +1 since first element is associated with state (i).
          termInd1 = termInd0+term{:}(1)-1;
          for ctrProd = 1 : 1 : numProd

            termSub = {subSum{ctrProd} subSTT{termInd0:termInd1}};
            prodTerm = STT{term{:}(ctrProd)}(termSub{:});
            prod = prod*prodTerm;
            
            if ctrProd ~= numProd
              termInd0 = termInd1+1;
              termInd1 = termInd0+term{:}(ctrProd+1)-1;
            end

          end

          % Multiply by factor to match coefficients in Taylor Series
          prod = prod*taylorSet{p,pJac}.multiplier{ctrTerm};

          % Add to terms within summation
          sumTerms = sumTerms + prod;

        end

        % Add term to summation
        indJac = {subSTT{1} subSum{:}};
        sumVal = sumVal + valJac{pJac}(indJac{:})*sumTerms;
        
      end

    end
    
    % Compute rate of change in STT term
    STTdot{p}(subSTT{:}) = sumVal;
  
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

