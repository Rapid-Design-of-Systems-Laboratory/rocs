function [] = writeSTT(fid,fidSymSTT,in,oc)
    fprintf(fidSymSTT,['function [STT] = ',in.robust.funcSymSTT,'()\n\n']);
    
    % Obtain terms for STT
    taylorSet = eomTermsSTT(in.cont.orderSTT);
    
    % Compute Jacobian matrices
    tempJac = oc.fAug;
    strJac = cell(orderSTT,1);
    for ctr = 1 : 1 : orderSTT
        
        strJac{ctr} = sym([]);
        
        % Determine size of Jacobian values are being assigned to
        sizeJac = oc.numAugStates*ones(1,ctr+1);
        
        % Construct empty cell array to store subscript indices
        sub = cell(1,length(sizeJac));
        
        % Compute Jacobian
        for ind = 1 : 1 : oc.numAugStates^(ctr+1)
            
            % Construct cell array of subscript elements in Jacobian
            [sub{:}] = ind2sub(sizeJac,ind);
            
            % Compute Jacobian
            strJac{ctr}(sub{:}) = diff(tempJac(sub{1:end-1}),oc.xAug(sub{end}));
            
        end
        
        % Save Jacobian
        tempJac = strJac{ctr};
        
    end
    
    % Determine number of elements in STTs
    numSTTelements = 0;
    for ctrOrder = 1 : 1 : orderSTT
        numSTTelements = numSTTelements + oc.numAugStates^(ctrOrder+1);
    end
    
    % Write to file
    if oc.numArcs > 1
        fprintf(fid,'function [xDotOut] = derivFuncSTT(t,X,region,p,const,constraint,scale,VEH,tRef,xRef) %%#ok<INUSD,INUSL>\n');
    else
        fprintf(fid,'function [xDotOut] = derivFuncSTT(t,X,p,const,constraint,scale,VEH,tRef,xRef) %%#ok<INUSD,INUSL>\n');
    end
    writeHeader(fid,in,out.VEH);
    fprintf(fid,'assert(isa(t, ''double''));\n');
    fprintf(fid,'assert(all(size(t)== [1 1]));\n');
    fprintf(fid,'assert(isa(X, ''double''));\n');
    fprintf(fid,['assert(all(size(X)== [',int2str(numSTTelements),' 1]));\n']);
    fprintf(fid,'assert(isa(p, ''double''));\n');
    fprintf(fid,['assert(all(size(p)== [',int2str(oc.numP + oc.numNu0 + oc.numNuF),' 1]));\n']);
    fprintf(fid,'assert(isa(tRef, ''double''));\n');
    fprintf(fid,'coder.varsize(''tRef'',[10000 1]);\n');
    fprintf(fid,'assert(isa(xRef, ''double''));\n');
    fprintf(fid,'coder.varsize(''xRef'',[10000 8]);\n');
    if oc.numArcs > 1
        fprintf(fid,'assert(isa(region, ''double''));\n');
    else
        fprintf(fid,'region = 1;\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'%%%% Parameters %%%%\n');
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'\n');
    % Time points
    for ctr = 1 : 1 : oc.numArcs
        fprintf(fid,['t',int2str(ctr),' = p(',int2str(ctr),');\n']);
    end
    % Costate jumps
    for ctr = 1 : 1 : length(oc.names)
        fprintf(fid,['pii_',oc.names{ctr},' = p(',int2str(2*ctr-1+oc.numArcs),':', ...
            int2str(2*ctr+oc.numArcs),');\n']);
    end
    pIndex = 2*length(oc.names)+oc.numArcs;
    fprintf(fid,'\n');
    % Initial constraint multipliers
    for ctr = 1 : 1 : oc.numNu0
        fprintf(fid,['nu0',int2str(ctr),' = p(',int2str(pIndex+ctr),');\n']);
    end
    fprintf(fid,'\n');
    % Terminal constraint multipliers
    for ctr = 1 : 1 : oc.numNuF
        fprintf(fid,['nuF',int2str(ctr),' = p(',int2str(pIndex+oc.numNu0+ctr),');\n']);
    end
    fprintf(fid,'\n');
    
    % Compute states and costates which are interpolated from reference
    % trajectory.
    fprintf(fid,'%% Compute states and costates at given time along reference trajectory\n');
    fprintf(fid,'xStar = interp1(tRef,xRef,t);\n');
    fprintf(fid,'%% States\n');
    for ctr = 1 : 1 : oc.numStates
        fprintf(fid,[char(oc.x(ctr)),' = xStar(',int2str(ctr),');\n']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'%% Costates\n');
    for ctr = 1 : 1 : oc.numStates
        fprintf(fid,[char(oc.lambda(ctr)),' = xStar(',int2str(ctr+oc.numStates),');\n']);
    end
    fprintf(fid,'\n');
    writeConstants(fid,in.const,false);
    fprintf(fid,'\n');
    writeScalingParameters(fid,in.scale,false);
    fprintf(fid,'\n');
    writeVehicleParameters(fid,out.VEH,false);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    % Assign costate derivatives based on trajectory arc
    fprintf(fid,'switch region\n');
    
    % Unconstrained arc
    for ctr = 1 : 1 : length(in.oc.sequence.unconstrained)
        arcNum = in.oc.sequence.unconstrained(ctr);
        fprintf(fid,['  case {',int2str(arcNum),'} %% unconstrained arc\n']);
        fprintf(fid,'     if lamGAM*Cl > 0\n');
        fprintf(fid,'       bank = pi;\n');
        fprintf(fid,'     else\n');
        fprintf(fid,'      bank = 0;\n');
        fprintf(fid,'     end\n');
        fprintf(fid,'    muu = 0;\n');
        muu = sym(0);
        if arcNum == 1
            timeShift = sym('t1');
        else
            timeShift = sym(['t',int2str(arcNum)]) - sym(['t',int2str(arcNum-1)]);
        end
        
        % Write equations to compute Jacobian
        startInd = 1;
        for p = 1 : 1 : orderSTT
            
            % Determine size of Jacobian values are being assigned to
            sizeJac = oc.numAugStates*ones(1,p+1);
            
            % Initialize Jacobian parameter for codegen
            strSize = '';
            for strCtr = 1 : 1 : p+1
                strSize = [strSize,int2str(oc.numAugStates)];
                if strCtr ~= p+1
                    strSize = [strSize,','];
                end
            end
            fprintf(fid,['valJac',int2str(p),' = nan(',strSize,');\n']);
            
            % Construct empty cell array to store subscript indices
            sub = cell(1,length(sizeJac));
            
            % Compute Jacobian
            for ind = 1 : 1 : oc.numAugStates^(p+1)
                
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
                fprintf(fid,['valJac',int2str(p),'(',strInd,') = ', ...
                    char(strJac{p}(sub{:})),';\n']);
                
                % Also write function to symbolically evaluate Jacobian
                fprintf(fidSymSTT,['valJac',int2str(p),'{',strInd,'} = sym(''', ...
                    char(strJac{p}(sub{:})),''');\n']);
                % Write symbolic STT elements
                fprintf(fidSymSTT,['STT',int2str(p),'{',strInd,'} = sym(''STT', ...
                    strSTT,''');\n']);
                
            end
            
            % Construct STT into matrix form for easy indexing
            dimArg = num2cell(sizeJac);
            fprintf(fid,['STT',int2str(p),' = reshape(X(',int2str(startInd),':', ...
                int2str(oc.numAugStates^(p+1)+startInd-1),',1),[',int2str([dimArg{:}]), ...
                ']);\n']);
            
            % Increment start index
            startInd = oc.numAugStates^(p+1)+startInd;
            
        end
        
        % Write equations of motion for STT terms
        
        % Compute STT rates for each order in Taylor Series approximation
        for p = 1 : 1 : orderSTT
            
            % Determine size and subscripts of STT
            sizeSTT = oc.numAugStates*ones(1,p+1);
            subSTT = cell(1,length(sizeSTT));
            
            % Initialize STTdot for codegen
            strSize = '';
            for strCtr = 1 : 1 : p+1
                strSize = [strSize,int2str(oc.numAugStates)];
                if strCtr ~= p+1
                    strSize = [strSize,','];
                end
            end
            fprintf(fid,['STTdot',int2str(p),' = nan(',strSize,');\n']);
            fprintf(fidSymSTT,['STT.STTdot',int2str(p),' = sym(nan(',strSize,'));\n']);
            
            % Write loops to go through each term in STT (a, b, c, etc)
            for loopCtr = 1 : 1 : p+1
                fprintf(fid,['for ctr',int2str(loopCtr),' = 1 : 1 : ', ...
                    int2str(oc.numAugStates),'\n']);
                fprintf(fidSymSTT,['for ctr',int2str(loopCtr),' = 1 : 1 : ', ...
                    int2str(oc.numAugStates),'\n']);
            end
            
            % Loop through each Jacobian term. Perform summation for term associated
            % with each Jacobian.
            fprintf(fid,'  sumVal = 0;\n');
            fprintf(fidSymSTT,'  sumVal = sym(''0'');\n');
            for pJac = 1 : 1 : p
                
                % Determine number of terms asocciated with Jacobian
                numTerms = length(taylorSet{p,pJac}.term);
                
                % Summation over alpha, beta, etc.
                fprintf(fid,['for sumCtr',int2str(pJac),' = 1 : 1 : ', ...
                    int2str(oc.numAugStates),'\n']);
                fprintf(fidSymSTT,['for sumCtr',int2str(pJac),' = 1 : 1 : ', ...
                    int2str(oc.numAugStates),'\n']);
                
                % Add terms within summation
                fprintf(fid,'sumVal = sumVal + (');
                fprintf(fidSymSTT,'sumVal = sumVal + (');
                for ctrTerm = 1 : 1 : numTerms
                    
                    % Get order for each product in term
                    term = taylorSet{p,pJac}.term(ctrTerm);
                    
                    % Determine number of products in each term
                    oc.numProd = length(term{:});
                    
                    % Perform multiplication of products
                    termInd0 = 2; % Shift by +1 since first element is associated with state (i).
                    termInd1 = termInd0+term{:}(1)-1;
                    for ctrProd = 1 : 1 : oc.numProd
                        
                        strInd = ['sumCtr',int2str(ctrProd)];
                        for strCtr = termInd0 : 1 : termInd1
                            strInd = [strInd,',ctr',int2str(strCtr)];
                        end
                        fprintf(fid,['STT',int2str(term{:}(ctrProd)),'(',strInd,')']);
                        fprintf(fidSymSTT,['STT',int2str(term{:}(ctrProd)),'{',strInd,'}']);
                        
                        if ctrProd ~= oc.numProd
                            termInd0 = termInd1+1;
                            termInd1 = termInd0+term{:}(ctrProd+1)-1;
                            fprintf(fid,'*');
                            fprintf(fidSymSTT,'*');
                        end
                        
                    end
                    
                    % Multiply by factor to match coefficients in Taylor Series
                    fprintf(fid,['*',num2str(taylorSet{p,pJac}.multiplier{ctrTerm})]);
                    fprintf(fidSymSTT,['*',num2str(taylorSet{p,pJac}.multiplier{ctrTerm})]);
                    
                    % Add to terms within summation
                    if ctrTerm ~= numTerms
                        fprintf(fid,'+');
                        fprintf(fidSymSTT,'+');
                    end
                    
                end
                
                fprintf(fid,')');
                fprintf(fidSymSTT,')');
                
                % Add term to summation
                strInd = 'ctr1';
                for ctrSum = 1 : 1 : pJac
                    strInd = [strInd,',sumCtr',int2str(ctrSum)];
                end
                fprintf(fid,['*valJac',int2str(pJac),'(',strInd,')']);
                fprintf(fidSymSTT,['*valJac',int2str(pJac),'{',strInd,'}']);
                
                fprintf(fid,';\n'); % End equation
                fprintf(fidSymSTT,';\n'); % End equation
                
            end
            
            % Close loops associated with summation (alpha, beta)
            for ctrSum = 1 : 1 : pJac
                fprintf(fid,'end\n');
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
            fprintf(fid,['STTdot',int2str(p),'(',strInd,') = real(sumVal);\n']);
            fprintf(fidSymSTT,['STT.STTdot',int2str(p),'(',strInd,') = sumVal;\n']);
            
            % Close loops associated with each element in STTs
            for loopCtr = 1 : 1 : p+1
                fprintf(fid,'end\n');
                fprintf(fidSymSTT,'end\n');
            end
            
        end
        
        % Reshape STT rates to column vector for output
        % Initialize Jacobian parameter for codegen
        numTerms = 0;
        for strCtr = 1 : 1 : p
            numTerms = numTerms + oc.numAugStates^(strCtr+1);
        end
        fprintf(fid,['xDotOut = nan(',int2str(numTerms),',1);\n']);
        
        startInd = 1;
        for ctrXdot = 1 : 1 : orderSTT
            numRows = oc.numAugStates^(ctrXdot+1);
            fprintf(fid,['xDotOut(',int2str(startInd),':', ...
                int2str(oc.numAugStates^(ctrXdot+1)+startInd-1),',1) = reshape(STTdot', ...
                int2str(ctrXdot),',',int2str(numRows),',1)*',char(timeShift),';\n']);
            startInd = startInd + oc.numAugStates^(ctrXdot+1);
        end
        
    end
    fprintf(fid,'  otherwise\n');
    fprintf(fid,'    error(''Incorrect region'');\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    fprintf(fid,'return\n');
    fprintf(fid,'\n');
    
    fprintf(fidSymSTT,'\n');
    fprintf(fidSymSTT,'return\n');
    fprintf(fidSymSTT,'\n');
end