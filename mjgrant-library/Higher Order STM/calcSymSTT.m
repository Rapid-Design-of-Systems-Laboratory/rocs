function [X0] = calcSymSTT(x,xDot,order,autocodeDir)

  % Open function that will be used to compute STTs symbolically
  fidSymSTT = fopen(['getSymSTT.m'],'w');
  fprintf(fidSymSTT,['function [STT,X0] = getSymSTT()\n\n']);

% Augmented state vector
numStates = length(x);

% Construct initial values for STT regardless if function is made.
for ctr = 1 : 1 : order

  % Determine size of Jacobian values are being assigned to
  sizeJac = numStates*ones(1,ctr+1);

  % Construct empty cell array to store subscript indices
  sub = cell(1,length(sizeJac));

  % Compute Jacobian
  for ind = 1 : 1 : numStates^(ctr+1)

    % Construct cell array of subscript elements in Jacobian    
    [sub{:}] = ind2sub(sizeJac,ind);

    % Construct initial values for STT. Identity for first order STT and zero 
    % for rest.
    if (ctr == 1) && (sum(diff([sub{:}])) == 0)
      val = 1;
    else
      val = 0;
    end
    X0{ctr}(ind,1) = val;
    
  end
  
end
  
% Obtain terms for STT
taylorSet = eomTermsSTT(order);
  
% Compute Jacobian matrices
tempJac = xDot;
strJac = cell(order,1);
for ctr = 1 : 1 : order
  
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
    
  end
  
  % Save Jacobian
  tempJac = strJac{ctr};
  
end
  
% Determine number of elements in STTs
numSTTelements = 0;
for ctrOrder = 1 : 1 : order
  numSTTelements = numSTTelements + numStates^(ctrOrder+1);
end


    
    % Write equations to compute Jacobian
    startInd = 1;
    for p = 1 : 1 : order

      % Determine size of Jacobian values are being assigned to
      sizeJac = numStates*ones(1,p+1);

      % Initialize Jacobian parameter for codegen
      strSize = '';
      for strCtr = 1 : 1 : p+1
        strSize = [strSize,int2str(numStates)];
        if strCtr ~= p+1
          strSize = [strSize,','];
        end
      end

      % Construct empty cell array to store subscript indices
      sub = cell(1,length(sizeJac));

      % Compute Jacobian
      for ind = 1 : 1 : numStates^(p+1)

        % Construct cell array of subscript elements in Jacobian    
        [sub{:}] = ind2sub(sizeJac,ind);

        % Write Jacobian
        strInd = '';
        strSTT = '';
        for strCtr = 1 : 1 : length([sub{:}])
          strInd = [strInd,int2str(sub{strCtr})];
          strSTT = [strSTT,int2str(sub{strCtr})];
          if strCtr ~= length([sub{:}])
            strInd = [strInd,','];
          end
        end
        
        % Also write function to symbolically evaluate Jacobian
        fprintf(fidSymSTT,['valJac',int2str(p),'{',strInd,'} = sym(''', ...
          char(strJac{p}(sub{:})),''');\n']);
        % Write symbolic STT elements
        fprintf(fidSymSTT,['STT',int2str(p),'{',strInd,'} = sym(''STT', ...
          strSTT,''');\n']);
        
      end
      
      % Construct STT into matrix form for easy indexing
      dimArg = num2cell(sizeJac);
      
      % Increment start index
      startInd = numStates^(p+1)+startInd;
      
    end
    
    % Write equations of motion for STT terms
    
    % Compute STT rates for each order in Taylor Series approximation
    for p = 1 : 1 : order
      
      % Determine size and subscripts of STT
      sizeSTT = numStates*ones(1,p+1);
      subSTT = cell(1,length(sizeSTT));
      
      % Initialize STTdot for codegen
      strSize = '';
      for strCtr = 1 : 1 : p+1
        strSize = [strSize,int2str(numStates)];
        if strCtr ~= p+1
          strSize = [strSize,','];
        end
      end
      fprintf(fidSymSTT,['STT.STTdot',int2str(p),' = sym(nan(',strSize,'));\n']);
      
      % Write loops to go through each term in STT (a, b, c, etc)
      for loopCtr = 1 : 1 : p+1
        fprintf(fidSymSTT,['for ctr',int2str(loopCtr),' = 1 : 1 : ', ...
          int2str(numStates),'\n']);
      end
      
      % Loop through each Jacobian term. Perform summation for term associated 
      % with each Jacobian.
      fprintf(fidSymSTT,'  sumVal = sym(''0'');\n');
      for pJac = 1 : 1 : p
        
        % Determine number of terms asocciated with Jacobian
        numTerms = length(taylorSet{p,pJac}.term);
        
        % Summation over alpha, beta, etc.
        fprintf(fidSymSTT,['for sumCtr',int2str(pJac),' = 1 : 1 : ', ...
          int2str(numStates),'\n']);
        
        % Add terms within summation
        fprintf(fidSymSTT,'sumVal = sumVal + (');
        for ctrTerm = 1 : 1 : numTerms
          
          % Get order for each product in term
          term = taylorSet{p,pJac}.term(ctrTerm);
          
          % Determine number of products in each term
          numProd = length(term{:});
          
          % Perform multiplication of products
          termInd0 = 2; % Shift by +1 since first element is associated with state (i).
          termInd1 = termInd0+term{:}(1)-1;
          for ctrProd = 1 : 1 : numProd
            
            strInd = ['sumCtr',int2str(ctrProd)];
            for strCtr = termInd0 : 1 : termInd1
              strInd = [strInd,',ctr',int2str(strCtr)];
            end
            fprintf(fidSymSTT,['STT',int2str(term{:}(ctrProd)),'{',strInd,'}']);
            
            if ctrProd ~= numProd
              termInd0 = termInd1+1;
              termInd1 = termInd0+term{:}(ctrProd+1)-1;
              fprintf(fidSymSTT,'*');
            end
            
          end
          
          % Multiply by factor to match coefficients in Taylor Series
          fprintf(fidSymSTT,['*',num2str(taylorSet{p,pJac}.multiplier{ctrTerm})]);
          
          % Add to terms within summation
          if ctrTerm ~= numTerms
            fprintf(fidSymSTT,'+'); 
          end
          
        end
        
        fprintf(fidSymSTT,')');
        
        % Add term to summation
        strInd = 'ctr1';
        for ctrSum = 1 : 1 : pJac
          strInd = [strInd,',sumCtr',int2str(ctrSum)];
        end
        fprintf(fidSymSTT,['*valJac',int2str(pJac),'{',strInd,'}']);
        
        fprintf(fidSymSTT,';\n'); % End equation
        
      end
      
      % Close loops associated with summation (alpha, beta)
      for ctrSum = 1 : 1 : pJac
        fprintf(fidSymSTT,'end\n');
      end
      
      % Assign summation to element in STT
      strInd = '';
      for loopCtr = 1 : 1 : p+1
        strInd = [strInd,'ctr',int2str(loopCtr)];
        if loopCtr ~= p+1
          strInd = [strInd,','];
        end
      end
      fprintf(fidSymSTT,['STT.STTdot',int2str(p),'(',strInd,') = sumVal;\n']);
      
      % Close loops associated with each element in STTs
      for loopCtr = 1 : 1 : p+1
        fprintf(fidSymSTT,'end\n');
      end
      
    end
    
    % Reshape STT rates to column vector for output
    % Initialize Jacobian parameter for codegen
    numTerms = 0;
    for strCtr = 1 : 1 : p
      numTerms = numTerms + numStates^(strCtr+1);
    end

    startInd = 1;
    for ctrXdot = 1 : 1 : order
      numRows = numStates^(ctrXdot+1);
      startInd = startInd + numStates^(ctrXdot+1);
    end
  
  fprintf(fidSymSTT,'\n');
  fprintf(fidSymSTT,'return\n');
  fprintf(fidSymSTT,'\n');
  fclose(fidSymSTT);
  
  end