function [] = writeDerivFuncJac(fid,in,oc,veh)
if oc.numArcs > 1
    fprintf(fid,'function [xDot,xDotParam] = derivFuncJac(t,X,region,p,const,constraint,scale,VEH,x0,xf) %%#ok<INUSD,INUSL>\n');
else
    fprintf(fid,'function [xDot,xDotParam] = derivFuncJac(t,X,p,const,constraint,scale,VEH,x0,xf) %%#ok<INUSD,INUSL>\n');
end
writeHeader(fid,in,veh);
fprintf(fid,'assert(isa(t, ''double''));\n');
fprintf(fid,'assert(all(size(t)== [1 1]));\n');
fprintf(fid,'assert(isa(X, ''double''));\n');
fprintf(fid,['assert(all(size(X)== [',int2str(2*oc.numStates),' 1]));\n']);
fprintf(fid,'assert(isa(p, ''double''));\n');
fprintf(fid,['assert(all(size(p)== [',int2str(oc.numP + oc.numNu0 + oc.numNuF),' 1]));\n']);
fprintf(fid,'assert(isa(x0, ''double''));\n');
fprintf(fid,['assert(all(size(x0)== [',int2str(oc.numStates),' 1]));\n']);
fprintf(fid,'assert(isa(xf, ''double''));\n');
fprintf(fid,['assert(all(size(xf)== [',int2str(oc.numStates),' 1]));\n']);
fprintf(fid,['xDot = nan(',num2str(2*oc.numStates),',',num2str(2*oc.numStates),');\n']);
fprintf(fid,['xDotParam = nan(',num2str(2*oc.numStates),',1);\n']);
if oc.numArcs > 1
    fprintf(fid,'assert(isa(region, ''double''));\n');
else
    fprintf(fid,'region = 1;\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%%%% Parameters %%%%\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
% Time points
for ctr = 1 : 1 : oc.numArcs
    fprintf(fid,['t',int2str(ctr),' = p(',int2str(ctr),');\n']);
    jacP = [jacP sym(['t',int2str(ctr)])];
end
% Costate jumps
for ctr = 1 : 1 : length(oc.names)
    fprintf(fid,['pii_',oc.names{ctr},' = p(',int2str(2*ctr-1+oc.numArcs),':', ...
        int2str(2*ctr+oc.numArcs),');\n']);
    jacP = [jacP sym(['pii_',oc.names{ctr},'1'])];
    jacP = [jacP sym(['pii_',oc.names{ctr},'2'])];
end
pIndex = 2*length(oc.names)+oc.numArcs;
fprintf(fid,'\n');
% Initial constraint multipliers
for ctr = 1 : 1 : oc.numNu0
    fprintf(fid,['nu0',int2str(ctr),' = p(',int2str(pIndex+ctr),');\n']);
    jacP = [jacP sym(['nu0',int2str(ctr)])];
end
fprintf(fid,'\n');
% Terminal constraint multipliers
for ctr = 1 : 1 : oc.numNuF
    fprintf(fid,['nuF',int2str(ctr),' = p(',int2str(pIndex+oc.numNu0+ctr),');\n']);
    jacP = [jacP sym(['nuF',int2str(ctr)])];
end
fprintf(fid,'\n');
fprintf(fid,'%% States\n');
for ctr = 1 : 1 : oc.numStates
    fprintf(fid,[char(oc.x(ctr)),' = X(',int2str(ctr),');\n']);
end
fprintf(fid,'\n');
fprintf(fid,'%% Costates\n');
for ctr = 1 : 1 : oc.numStates
    fprintf(fid,[char(oc.lambda(ctr)),' = X(',int2str(ctr+oc.numStates),');\n']);
end
fprintf(fid,'\n');
writeConstants(fid,in.const,false);
fprintf(fid,'\n');
writeScalingParameters(fid,in.scale,false);
fprintf(fid,'\n');
writeVehicleParameters(fid,veh,false);
fprintf(fid,'\n');

% Assign costate derivatives based on trajectory arc
fprintf(fid,'switch region\n');
% Unconstrained arc
for ctr = 1 : 1 : length(in.oc.sequence.unconstrained)
    arcNum = in.oc.sequence.unconstrained(ctr);
    fprintf(fid,['\n  case {',int2str(arcNum),'} %% unconstrained arc\n']);
    if in.oc.sequence.singularArc(ctr)
        
        fprintf(fid,['     quantity = ',char(bankCosArg),';\n']);
        fprintf(fid,'     if abs(quantity) > 1\n');
        fprintf(fid,'       quantity = sign(quantity);\n');
        fprintf(fid,'     end\n');
        fprintf(fid,'     bank = acos(quantity);\n');
        
    else
        
        if in.oc.switch.enforced(ctr == in.oc.sequence.unconstrained)
            fprintf(fid,'     if lamGAM*Cl > 0\n');
            fprintf(fid,'       bank = pi;\n');
            fprintf(fid,'     else\n');
            fprintf(fid,'      bank = 0;\n');
            fprintf(fid,'     end\n');
        else
            if in.oc.switch.direction(ctr == in.oc.sequence.unconstrained) > 0
                fprintf(fid,'     bank = 0;\n');
            elseif in.oc.switch.direction(ctr == in.oc.sequence.unconstrained) < 0
                fprintf(fid,'     bank = pi;\n');
            else
                error('Switch direction not a correct value.');
            end
        end
        
    end
    fprintf(fid,'    muu = 0;\n');
    muu = sym(0);
    if arcNum == 1
        timeShift = sym('t1');
    else
        timeShift = sym(['t',int2str(arcNum)]) - sym(['t',int2str(arcNum-1)]);
    end
    varString = ['    xDot(',int2str(oc.numStates+1),':',int2str(2*oc.numStates),',:)'];
    func = oc.lambdaDot.unconstrained*timeShift;
    writeJacobian(fid,varString,func,jacY);
    fprintf(fid,'\n');
    % Jacobian with respect to parameters
    varString = ['    xDotParam(',int2str(oc.numStates+1),':',int2str(2*oc.numStates),',:)'];
    writeJacobian(fid,varString,func,jacP);
    fprintf(fid,'\n');
    fprintf(fid,'%% Equations of motion\n');
    varString = ['xDot(1:',int2str(oc.numStates),',:)'];
    func = oc.xDot*timeShift;
    writeJacobian(fid,varString,func,jacY);
    varString = ['xDotParam(1:',int2str(oc.numStates),',:)'];
    writeJacobian(fid,varString,func,jacP);
    
end

% Constrained arcs
for ctr1 = 1 : 1 : length(oc.names)
    % Write case statement for constraint
    for ctr2 = 1 : 1 : length(in.oc.sequence.(oc.names{ctr1}))
        arcNum = in.oc.sequence.(oc.names{ctr1})(ctr2);
        fprintf(fid,['\n  case {',int2str(arcNum),'} %% ',oc.names{ctr1},' arc\n']);
        fprintf(fid,'Cl = ClConstraint;\n');
        % Write control values
        for ctr3 = 1 : 1 : length(in.oc.u)
            fprintf(fid,['    ',char(in.oc.u(ctr3)),' = ',char(oc.uOptimal.(oc.names{ctr1})(ctr3)),';\n']);
        end
        % Write muu values
        muConstraint = sym(['mu_',oc.names{ctr1}]);
        fprintf(fid,['    ',char(muConstraint),' = ',char(oc.muu.(oc.names{ctr1})),';\n']);
        if arcNum == 1
            timeShift = sym('t1');
        else
            timeShift = sym(['t',int2str(arcNum)]) - sym(['t',int2str(arcNum-1)]);
        end
        varString = ['    xDot(',int2str(length(oc.x)+1),':',int2str(2*length(oc.x)),',:)'];
        func = subs(oc.lambdaDot.(oc.names{ctr1})(:)*timeShift,{in.oc.u(1),muConstraint},{oc.uOptimal.(oc.names{ctr1})(1),oc.muu.(oc.names{ctr1})});
        writeJacobian(fid,varString,func,jacY);
        fprintf(fid,'\n');
        % Write derivative with respect to parameters
        varString = ['    xDotParam(',int2str(oc.numStates+1),':',int2str(2*oc.numStates),',:)'];
        writeJacobian(fid,varString,func,jacP);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        
        fprintf(fid,'%% Equations of motion\n');
        varString = ['xDot(1:',int2str(oc.numStates),',:)'];
        func = subs(oc.xDot*timeShift,{in.oc.u(1),muConstraint},{oc.uOptimal.(oc.names{ctr1})(1),oc.muu.(oc.names{ctr1})});
        writeJacobian(fid,varString,func,jacY);
        varString = ['xDotParam(1:',int2str(oc.numStates),',:)'];
        writeJacobian(fid,varString,func,jacP);
        
    end
end
fprintf(fid,'  otherwise\n');
fprintf(fid,'    error(''Incorrect region'');\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');
fprintf(fid,'return\n');
fprintf(fid,'\n');