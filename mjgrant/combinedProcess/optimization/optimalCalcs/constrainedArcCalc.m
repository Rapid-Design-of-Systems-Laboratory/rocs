function [oc2] = constrainedArcCalc(in,oc,veh)

	% Loop through each constraint
	for ctr1 = 1 : 1 : length(oc.names)
	    % Take derivatives until control found
	    index = 1;
	    while true
        
	        if isempty(strfind(char(oc.S.(oc.names{ctr1})(index,1)),char(in.oc.u)))
	            % Control not found yet, take another derivative
	            oc.S.(oc.names{ctr1})(index+1,1) = jacobian(oc.S.(oc.names{ctr1})(index,1),oc.x)*oc.xDot;
	            index = index + 1; % Update the index
	        else
	            % Control found, exit
	            break;
	        end
        
	    end
	    % Get dN/dx. Includes all Jacobians up to but not including set with control
	    % variable present.
	    for ctr2 = 1 : 1 : index - 1
	        oc.dNdx.(oc.names{ctr1})(1:length(oc.x),ctr2) = jacobian(oc.S.(oc.names{ctr1})(ctr2,1),oc.x).';
	    end
    
	    % Write external function for dNdx
	    funcName = [oc.names{ctr1},'_dNdx'];
	    fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
	    % fprintf(fid,['function dNdx = ',funcName,'(r,v,gam,const,constraint,scale,VEH)\n\n']);
	    fprintf(fid,['function dNdx = ',funcName,'(']);
	    % Add states as function inputs
	    for row = 1 : 1 : oc.numStates
	        fprintf(fid,'%s',char(oc.x(row,1)));
	        if row ~= oc.numStates
	            fprintf(fid,',');
	        end
	    end
	    % Add rest of the parameters
	    fprintf(fid,',const,constraint,scale,VEH)\n\n');
	    writeHeader(fid,in,veh);
	    fprintf(fid,'\n');
	    writeConstants(fid,in.const,false);
	    fprintf(fid,'\n');
	    % Generalize this (remove alfa)
	    fprintf(fid,'syms alfa\n');
	    writeVehicleParameters(fid,veh,false);
    
	    fprintf(fid,'\n');
	    [numRow,numCol] = size(oc.dNdx.(oc.names{ctr1}));
	    fprintf(fid,'dNdx = [');
	    for row = 1 : 1 : numRow
	        for col = 1 : 1 : numCol
	            fprintf(fid,'%s',char(oc.dNdx.(oc.names{ctr1})(row,col)));
	            if col ~= numCol
	                fprintf(fid,', ');
	            elseif row ~= numRow
	                fprintf(fid,';\n  ');
	            end
	        end
	    end
	    fprintf(fid,'];\n\n');
	    fprintf(fid,'return\n\n');
	    fclose(fid);
    
	    % Get dN/dt
	    for ctr2 = 1 : 1 : index - 1
	        oc.dNdt.(oc.names{ctr1})(ctr2,1) = jacobian(oc.S.(oc.names{ctr1})(ctr2,1),sym('t'));
	    end
    
	    % Write external function for dNdt
	    funcName = [oc.names{ctr1},'_dNdt'];
	    fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
	    % fprintf(fid,['function dNdt = ',funcName,'(r,v,gam,const,constraint,scale,VEH)\n\n']);
	    fprintf(fid,['function dNdt = ',funcName,'(']);
	    % Add states as function inputs
	    for row = 1 : 1 : oc.numStates
	        fprintf(fid,'%s',char(oc.x(row,1)));
	        if row ~= oc.numStates
	            fprintf(fid,',');
	        end
	    end
	    % Add rest of the parameters
	    fprintf(fid,',const,constraint,scale,VEH)\n\n');
	    writeHeader(fid,in,veh);
	    fprintf(fid,'\n');
	    writeConstants(fid,in.const,false);
	    fprintf(fid,'\n');
	    % Generalize this (remove alfa)
	    fprintf(fid,'syms alfa\n');
	    writeVehicleParameters(fid,veh,false);
	    fprintf(fid,'\n');
	    [numRow,numCol] = size(oc.dNdt.(oc.names{ctr1}));
	    fprintf(fid,'dNdt = [');
	    for row = 1 : 1 : numRow
	        for col = 1 : 1 : numCol
	            fprintf(fid,'%s',char(oc.dNdt.(oc.names{ctr1})(row,col)));
	            if col ~= numCol
	                fprintf(fid,', ');
	            elseif row ~= numRow
	                fprintf(fid,';\n  ');
	            end
	        end
	    end
	    fprintf(fid,'];\n\n');
	    fprintf(fid,'return\n\n');
	    fclose(fid);
    
	    % Solve for bank along constraint
	    % Generalize this (remove bank)
	    bankMatrix = solve([char(oc.S.(oc.names{ctr1})(end)),'=0'],char(in.oc.u));
	    % Select proper index from equation. Need to generalize this.
	    if ctr1 == 1
	        matrixIndex = 1;
	    else
	        matrixIndex = 2;
	    end
	    oc.uOptimal.(oc.names{ctr1}) = bankMatrix(matrixIndex);
    
	    % Write external function for optimal bank along constraint
	    funcName = [oc.names{ctr1},'_control'];
	    fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
	    % fprintf(fid,['function control = ',funcName,'(r,v,gam,const,constraint,scale,VEH)\n\n']);
	    fprintf(fid,['function control = ',funcName,'(']);
	    % Add states as function inputs
	    for row = 1 : 1 : oc.numStates
	        fprintf(fid,'%s',char(oc.x(row,1)));
	        if row ~= oc.numStates
	            fprintf(fid,',');
	        end
	    end
	    % Add rest of the parameters
	    fprintf(fid,',const,constraint,scale,VEH)\n\n');
	    fprintf(fid,'\n');
	    writeConstants(fid,in.const,false);
	    fprintf(fid,'\n');
	    % Generalize this (remove alfa)
	    fprintf(fid,'syms alfa\n');
	    writeVehicleParameters(fid,veh,false);
	    fprintf(fid,'\n');
	    [numRow,numCol] = size(oc.uOptimal.(oc.names{ctr1}));
	    fprintf(fid,'dNdt = [');
	    for row = 1 : 1 : numRow
	        for col = 1 : 1 : numCol
	            fprintf(fid,'%s',char(oc.uOptimal.(oc.names{ctr1})(row,col)));
	            if col ~= numCol
	                fprintf(fid,', ');
	            elseif row ~= numRow
	                fprintf(fid,';\n  ');
	            end
	        end
	    end
	    fprintf(fid,'];\n\n');
	    fprintf(fid,'return\n\n');
	    fclose(fid);
    
	    % Compute Hamiltonian along constraint
	    muConstraint = sym(['mu_',oc.names{ctr1}]);
	    oc.H.(oc.names{ctr1}) = in.oc.Lagrange + oc.lambda.'*oc.xDot + muConstraint*oc.S.(oc.names{ctr1})(end,1);
    
	    % Compute partial of Hamiltonian with respect to control
	    oc.Hu.(oc.names{ctr1}) = diff(oc.H.(oc.names{ctr1}),in.oc.u);
    
	    % Compute Lagrange multiplier derivatives
	    oc.lambdaDot.(oc.names{ctr1}) = -jacobian(oc.H.(oc.names{ctr1}),oc.x).';
    
	    % Calculate constraint Lagrange multiplier
	    oc.muu.(oc.names{ctr1}) = solve([char(oc.Hu.(oc.names{ctr1})),'=0'],char(muConstraint));
	
	end
	oc2 = oc;
end