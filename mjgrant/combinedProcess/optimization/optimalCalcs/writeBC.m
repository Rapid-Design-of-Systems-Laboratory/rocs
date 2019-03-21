function [] = writeBC(fid,in,oc)

%%%%%%%%%%%%%%%%%%
%% Write Header %%
%%%%%%%%%%%%%%%%%%

% if in.convertParametersToStates
% 	fprintf(fid,'function [zeroVec] = bc(YL,YR,p,const,constraint,arcSequence,interiorPointConstraintSequence,interiorPointNumLagrangeMultipliers,x0,xf)\n');
% else
fprintf(fid,'function [zeroVec] = bc(YL,YR,p,const,constraint,arcSequence,interiorPointConstraintSequence,interiorPointNumLagrangeMultipliers,x0,xf)\n');
% end

% if in.rootSolving == 0
writeHeader(fid,in);

if ~in.convertParametersToStates
    fprintf(fid,['coder.varsize(''YL'', [',int2str(oc.num.states*2),' ',int2str(in.maxNumArcs),']);\n']);
    fprintf(fid,['coder.varsize(''YR'', [',int2str(oc.num.states*2),' ',int2str(in.maxNumArcs),']);\n']);
else
    fprintf(fid,['coder.varsize(''YL'', [',int2str(oc.num.states*2+in.maxNumArcs),' ',int2str(in.maxNumArcs),']);\n']);
    fprintf(fid,['coder.varsize(''YR'', [',int2str(oc.num.states*2+in.maxNumArcs),' ',int2str(in.maxNumArcs),']);\n']);
end
fprintf(fid,'assert(isa(YL,''double''));\n');
fprintf(fid,'assert(isa(YR,''double''));\n');
fprintf(fid,'assert(isa(p,''double''));\n');
fprintf(fid,['coder.varsize(''p'', [',int2str(in.maxNumArcs + in.maxNumArcs*oc.num.states + oc.num.constraints.initial + oc.num.constraints.terminal),' 1]);\n']);
% fprintf(fid,['coder.varsize(''p'');\n']);

fprintf(fid,'assert(isa(x0,''double''));\n');
fprintf(fid,['assert(all(size(x0)== [',int2str(oc.num.states),' 1]));\n']);
fprintf(fid,'assert(isa(xf,''double''));\n');
fprintf(fid,['assert(all(size(xf)== [',int2str(oc.num.states),' 1]));\n']);
fprintf(fid,'assert(isa(arcSequence, ''double''));\n');
fprintf(fid,['coder.varsize(''arcSequence'', [1 ',int2str(in.maxNumArcs),']);\n']);
% fprintf(fid,'assert(isa(piiIndex, ''double''));\n');
% fprintf(fid,['coder.varsize(''piiIndex'', [',int2str(in.maxNumArcs*oc.num.states),',2]);\n']);

fprintf(fid,'assert(isa(interiorPointConstraintSequence, ''double''));\n');

if oc.num.constraints.path+oc.num.constraints.interiorPoint ~= 0
    fprintf(fid,['coder.varsize(''interiorPointConstraintSequence'', [',int2str(oc.num.constraints.path+oc.num.constraints.interiorPoint),' ',int2str(in.maxNumArcs),']);\n']);
else
    fprintf(fid,['assert(all(size(interiorPointConstraintSequence) == [0 0]));\n']);
end

fprintf(fid,'assert(isa(interiorPointNumLagrangeMultipliers, ''double''));\n');
if oc.num.lagrangeMultipliers.interiorPoint ~= 0
    fprintf(fid,['coder.varsize(''interiorPointNumLagrangeMultipliers'', [1 ',int2str(length(oc.num.lagrangeMultipliers.interiorPoint)),']);\n']);
else
    fprintf(fid,['assert(all(size(interiorPointNumLagrangeMultipliers) == [0 0]));\n']);
end


fprintf(fid,['coder.varsize(''constraintDeriv'', [',int2str(in.maxNumArcs*oc.num.states),' 1]);\n']);
fprintf(fid,['coder.varsize(''constraintDeriv1'', [',int2str(in.maxNumArcs*oc.num.states),' 1]);\n']);
fprintf(fid,['coder.varsize(''H_t_plus'', [',int2str(in.maxNumArcs),' 1]);\n']);
fprintf(fid,['coder.varsize(''H_t_minus'', [',int2str(in.maxNumArcs),' 1]);\n']);
fprintf(fid,['coder.varsize(''Harcs'', [',int2str(in.maxNumArcs),' 1]);\n']);
fprintf(fid,['coder.varsize(''continuityStates'', [',int2str(in.maxNumArcs*oc.num.states),' 1]);\n']);
fprintf(fid,['coder.varsize(''costateArcs'', [',int2str(in.maxNumArcs*oc.num.states),' 1]);\n']);
fprintf(fid,['coder.varsize(''dNdt'', [',int2str(in.maxNumArcs*oc.num.states),',1]);\n']);
fprintf(fid,['coder.varsize(''dNdx'', [',int2str(in.maxNumArcs*oc.num.states),',',int2str(oc.num.states),']);\n']);

fprintf(fid,['coder.varsize(''hamiltonianDiscontinuity'', [',int2str(in.maxNumArcs*oc.num.states),',1]);\n']);
fprintf(fid,['coder.varsize(''costateDiscontinuity'', [',int2str(in.maxNumArcs*oc.num.states),',',int2str(oc.num.states),']);\n']);

if in.verbose
    fprintf(fid,'coder.extrinsic(''keyboard'');\n');
end

fprintf(fid,'\n');

%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Number of arcs
fprintf(fid,'numArcs = length(arcSequence);\n');

fprintf(fid,'interiorPointConstraintIndex = 0;\n');
fprintf(fid,'lmInteriorPointIndex = 1;\n');

fprintf(fid,'interiorPointConstraintIndex1 = 0;\n');
fprintf(fid,'lmInteriorPointIndex1 = 1;\n');

fprintf(fid,'indexConstraintDeriv = 0;\n'); % used to keep track of number of tangency constraints
fprintf(fid,'indexConstraintDeriv1 = 0;\n');

fprintf(fid,'indexContinuityStates = 1;\n'); % used to keep track of number of state continuity equations across arcs
fprintf(fid,'indexCostateArcs = 1;\n'); % used to keep track of costate conditions across arcs

fprintf(fid,'indexHarcs = 0;\n'); % used to keep track of Hamiltonian conditions across arcs
fprintf(fid,'H_t_plus = NaN(numArcs,1);\n');
fprintf(fid,'H_t_minus = NaN(numArcs,1);\n');

if in.maxNumArcs > 1
    if ~in.convertParametersToStates
        fprintf(fid,['continuityStates = NaN((numArcs-1)*',int2str(oc.num.states),',1);\n']);
    else
        fprintf(fid,['continuityStates = NaN((numArcs-1)*(',int2str(oc.num.states),'+numArcs),1);\n']);
    end
    fprintf(fid,['costateArcs = NaN((numArcs-1)*',int2str(oc.num.states),',1);\n']);
    fprintf(fid,'Harcs = NaN(numArcs-1,1);\n');
else
    fprintf(fid,'Harcs = NaN(1,1);\n');
end

for ctrState = 1 : 1 : oc.num.states
    fprintf(fid,[char(oc.state.var(ctrState)),' = NaN;\n']);
end

for ctrCostate = 1 : 1 : oc.num.states
    fprintf(fid,[char(oc.costate.var(ctrCostate)),' = NaN;\n']);
end

for ctrControl = 1 : 1 : length(oc.control.var)
    fprintf(fid,[char(oc.control.var(ctrControl)),' = NaN;\n']);
end

if isfield(in,'quantity')
    names = fieldnames(in.quantity);
else
    names = {};
end

for ctrQuantity = 1 : 1 : length(names)
    fprintf(fid,[char(names{ctrQuantity}),' = NaN;\n']);
end

fprintf(fid,['xAndLambda0Constraint = NaN(',int2str(oc.num.constraints.initial+oc.num.states),',1);\n']);
fprintf(fid,['xAndLambdaFConstraint = NaN(',int2str(oc.num.constraints.terminal+oc.num.states),',1);\n']);

if ~in.convertParametersToStates
    fprintf(fid,['xAndLambda = NaN(',int2str(2*oc.num.states),',1);\n']);
else
    fprintf(fid,['xAndLambda = NaN(',int2str(2*oc.num.states+in.maxNumArcs),',1);\n']);
end
fprintf(fid,'pii = NaN;\n');
fprintf(fid,'pii1 = NaN;\n');

%%%%%%%%%%%%%%%%%%%%%%
%% Write Parameters %%
%%%%%%%%%%%%%%%%%%%%%%



fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%%%% Parameters %%%%\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');

% fprintf(fid,['p\n']);

% Time points
if ~in.convertParametersToStates
    fprintf(fid,['tSet = p(1:numArcs);\n']);
else
    %		fprintf(fid,['t_ = X(',num2str(oc.num.origStates),'+region,1);']);
end
fprintf(fid,'\n');

% Initial constraint multipliers
if ~in.convertParametersToStates
    for ctrNu0 = 1 : 1 : oc.num.constraints.initial
        fprintf(fid,['lmInitial',int2str(ctrNu0),' = p(numArcs+',int2str(ctrNu0),');\n']);
    end
else
    for ctrNu0 = 1 : 1 : oc.num.constraints.initial
        fprintf(fid,['lmInitial',int2str(ctrNu0),' = p(0+',int2str(ctrNu0),');\n']);
    end
end
fprintf(fid,'\n');
% Terminal constraint multipliers
if ~in.convertParametersToStates
    for ctrNuF = 1 : 1 : oc.num.constraints.terminal
        fprintf(fid,['lmTerminal',int2str(ctrNuF),' = p(numArcs+',int2str(oc.num.constraints.initial),'+',int2str(ctrNuF),');\n']);
    end
else
    for ctrNuF = 1 : 1 : oc.num.constraints.terminal
        fprintf(fid,['lmTerminal',int2str(ctrNuF),' = p(0+',int2str(oc.num.constraints.initial),'+',int2str(ctrNuF),');\n']);
    end
end
fprintf(fid,'\n');

% Costate jumps
if ~in.convertParametersToStates 
% fprintf(fid,'if numArcs >1\n');
     fprintf(fid,['lmInteriorPoint = p(numArcs+',int2str(oc.num.constraints.initial),'+',int2str(oc.num.constraints.terminal),'+1:end);\n']);
fprintf(fid,['lmInteriorPoint1 = 0;\n']);
% fprintf(fid,'else\n');
% fprintf(fid,['lmInteriorPoint = p(numArcs+',int2str(oc.num.constraints.initial),'+',int2str(oc.num.constraints.terminal),'+1:numArcs+',int2str(oc.num.constraints.initial),'+',int2str(oc.num.constraints.terminal),'+3)\n']);
% fprintf(fid,['lmInteriorPoint1 = p(numArcs+',int2str(oc.num.constraints.initial),'+',int2str(oc.num.constraints.terminal),'+3:end)\n']);
% fprintf(fid,'end\n');
else
    fprintf(fid,['lmInteriorPoint = p(0+',int2str(oc.num.constraints.initial),'+',int2str(oc.num.constraints.terminal),'+1:end)\n']);
end
fprintf(fid,'\n');
fprintf(fid,'dNdt = NaN(length(lmInteriorPoint),1);\n');
fprintf(fid,['dNdx = NaN(',int2str(oc.num.states),',length(lmInteriorPoint));\n']);
fprintf(fid,'constraintDeriv = NaN(length(lmInteriorPoint),1);\n');
fprintf(fid,'constraintDeriv1 = NaN(length(lmInteriorPoint1),1);\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write Constants and Constraint Values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(in,'const')
    writeConstants(fid,in.const,false);
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
if isfield(in,'constraintVal')
    writeConstraints(fid,in.constraintVal,false);
end
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Through Each Arc and Obtain Boundary Condition Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Also loop through left and right endpoint to reduce repeated code
fprintf(fid,'for ctrArc = 1 : 1 : numArcs\n');

fprintf(fid,'\tindexArc = arcSequence(ctrArc);\n');
fprintf(fid,'\thamiltonianDiscontinuity = zeros(1,1);\n');
fprintf(fid,['\tcostateDiscontinuity = zeros(',int2str(oc.num.states),',1);\n']);

fprintf(fid,'\thamiltonianDiscontinuity1 = zeros(1,1);\n');
fprintf(fid,['\tcostateDiscontinuity1 = zeros(',int2str(oc.num.states),',1);\n']);

fprintf(fid,'\tfor ctrEndpoint = 1 : 1 : 2\n');

fprintf(fid,'\t\tif ctrEndpoint == 1\n'); % Write states and costates from left boundary

fprintf(fid,'\t\t\t%% Left endpoint\n');
fprintf(fid,'\t\t\t%% States\n');
for ctrState = 1 : 1 : oc.num.states
    fprintf(fid,['\t\t\t',char(oc.state.var(ctrState)),' = YL(',int2str(ctrState),',ctrArc);\n']);
end
fprintf(fid,'\n');

fprintf(fid,'\t\t\t%% Costates\n');
for ctrCostate = 1 : 1 : oc.num.states
    fprintf(fid,['\t\t\t',char(oc.costate.var(ctrCostate)),' = YL(',int2str(ctrCostate+oc.num.states),',ctrArc);\n']);
end
fprintf(fid,'\n');

if in.convertParametersToStates
    fprintf(fid,'\t\t\t%% Independant Variables\n');
    fprintf(fid,'\t\t\tif(~isa(YL,''sym''))\n');
    fprintf(fid,['\t\t\t\ttSet = NaN(',int2str(in.maxNumArcs),',1);\n']);
    fprintf(fid,'\t\t\tend;\n');
    fprintf(fid,['\t\t\tfor ctr = 1 : 1 : numArcs\n']);
    fprintf(fid,['\t\t\t\ttSet(ctr) = YL(',int2str(2*oc.num.states),'+ctr,ctrArc);\n']);
    fprintf(fid,['\t\t\tend\n']);
    fprintf(fid,'\n');
end

% Write control input for boundary point
fprintf(fid,'\t\t\txAndLambda = [');
for ctrState = 1 : 1 : oc.num.states
    fprintf(fid,[char(oc.state.var(ctrState)),';']);
end
for ctrCostate = 1 : 1 : oc.num.states
    fprintf(fid,char(oc.costate.var(ctrCostate)));
    if ctrCostate ~= oc.num.states
        fprintf(fid,';');
    else
        if in.convertParametersToStates
            for ctr = 1 : 1 : in.maxNumArcs
                fprintf(fid,';NaN');
            end
        end
        fprintf(fid,'];\n');
    end
end
if in.convertParametersToStates
    fprintf(fid,'\t\t\tfor ctr = 1 : 1 : numArcs\n');
    fprintf(fid,['\t\t\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
    fprintf(fid,'\t\t\tend\n');
end

% Write derived quantities
if isfield(in,'quantity')
    writeQuantities(fid,in.quantity);
    fprintf(fid,'\n');
end

% Initial state and costate constraint
fprintf(fid,'\t\t\tif ctrArc == 1\n');
fprintf(fid,'\t\t\t\txAndLambda0Constraint = [\n');
for ctrState = 1 : 1 : oc.num.states
    % 	if ismember(oc.state.var(ctrState),oc.constraint.initial) % Initial state specified
    % 		fprintf(fid,['\t\t\t\t\t',char(oc.state.var(ctrState)),' - x0(',int2str(ctrState),');\n']);
    % 	end
    fprintf(fid,['\t\t\t\t\t',char(oc.costate.var(ctrState)),' - (',char(oc.costate.initial(ctrState)),');\n']);
end

for ctrInitialConstraint = 1 : 1 : length(oc.constraint.initial)
    fprintf(fid,['\t\t\t\t\t',char(oc.constraint.initial(ctrInitialConstraint)),';\n']);
end

fprintf(fid,'\t\t\t\t];\n');
fprintf(fid,'\t\t\tend\n');
fprintf(fid,'\t\telseif ctrEndpoint == 2\n'); % Write states and costates from right boundary

fprintf(fid,'\t\t\t%% Right endpoint\n');
fprintf(fid,'\t\t\t%% States\n');
for ctrStates = 1 : 1 : oc.num.states
    fprintf(fid,['\t\t\t',char(oc.state.var(ctrStates)),' = YR(',int2str(ctrStates),',ctrArc);\n']);
end

fprintf(fid,'\n');

fprintf(fid,'\t\t\t%% Costates\n');
for ctrCostates = 1 : 1 : oc.num.states
    fprintf(fid,['\t\t\t',char(oc.costate.var(ctrCostates)),' = YR(',int2str(ctrCostates+oc.num.states),',ctrArc);\n']);
end
fprintf(fid,'\n');

if in.convertParametersToStates
    fprintf(fid,'\t\t\t%% Independant Variables\n');
    fprintf(fid,'\t\t\tif(~isa(YR,''sym''))\n');
    fprintf(fid,['\t\t\t\ttSet = NaN(',int2str(in.maxNumArcs),',1);\n']);
    fprintf(fid,'\t\t\tend;\n');
    fprintf(fid,['\t\t\tfor ctr = 1 : 1 : numArcs\n']);
    fprintf(fid,['\t\t\t\ttSet(ctr) = YL(',int2str(2*oc.num.states),'+ctr,ctrArc);\n']);
    fprintf(fid,['\t\t\tend\n']);
    fprintf(fid,'\n');
end

% Write control input for boundary point
fprintf(fid,'\t\t\txAndLambda = [');
for ctrState = 1 : 1 : oc.num.states
    fprintf(fid,[char(oc.state.var(ctrState)),';']);
end
for ctrCostate = 1 : 1 : oc.num.states
    fprintf(fid,char(oc.costate.var(ctrCostate)));
    if ctrCostate ~= oc.num.states
        fprintf(fid,';');
    else
        if in.convertParametersToStates
            for ctr = 1 : 1 : in.maxNumArcs
                fprintf(fid,';NaN');
            end
        end
        fprintf(fid,'];\n');
    end
end
if in.convertParametersToStates
    fprintf(fid,'\t\t\tfor ctr = 1 : 1 : numArcs\n');
    fprintf(fid,['\t\t\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
    fprintf(fid,'\t\t\tend\n');
end

% Write derived quantities
if isfield(in,'quantity')
    writeQuantities(fid,in.quantity);
    fprintf(fid,'\n');
end

% Final state and costate constraint
fprintf(fid,'\t\t\tif ctrArc == numArcs\n');
fprintf(fid,'\t\t\t\txAndLambdaFConstraint = [\n');
for ctrState = 1 : 1 : oc.num.states
    % 	if ismember(oc.state.var(ctrState),oc.constraint.terminal) % Final state specified
    % 		fprintf(fid,['\t\t\t\t\t',char(oc.state.var(ctrState)),' - xf(',int2str(ctrState),');\n']);
    % 	end
    fprintf(fid,['\t\t\t\t\t',char(oc.costate.var(ctrState)),' - (',char(oc.costate.terminal(ctrState)),');\n']);
end
if isfield(in.oc.constraint,'terminal')
    for ctrTerminalConstraint = 1 : 1 : length(oc.constraint.terminal)
        fprintf(fid,['\t\t\t\t\t',char(oc.constraint.terminal(ctrTerminalConstraint)),';\n']);
    end
end
fprintf(fid,'\t\t\t\t];\n');
fprintf(fid,'\t\t\tend\n');
fprintf(fid,'\t\tend\n'); % if ctrEndpoint
fprintf(fid,'\n');

% Compute control
fprintf(fid,'\t\tswitch indexArc\n');

fprintf(fid,'\t\tcase {0} %% unconstrained arc\n');

% Create list of control names and symbols
controlNames = '';
controlSymNames = '';
for ctrControl = 1 : 1 : oc.num.controls
    controlNames = strcat(controlNames,char(oc.control.var(ctrControl)));
    controlSymNames = strcat(controlSymNames,['sym(''',char(oc.control.var(ctrControl)),''')']);
    if ctrControl ~= oc.num.controls
        controlNames = strcat(controlNames,',');
        controlSymNames = strcat(controlSymNames,',');
    end
end

% Use symbolic control to allow for external control function in GPU code
fprintf(fid,'\t\t\tif(~isa(YL,''sym''))\n');
fprintf(fid,'\t\t\t\t[');
fprintf(fid,controlNames);
fprintf(fid,',d2Hdu2] = computeControlUnconstrained(xAndLambda,const,constraint,numArcs);\n');
fprintf(fid,'\t\t\telse\n\t');
fprintf(fid,['\t\t\t[',controlNames,',d2Hdu2] = ']);
fprintf(fid,['deal(',controlSymNames,',0);\n']);
fprintf(fid,'\t\t\tend\n');
fprintf(fid,'\n');

for ctrConstraint = 1 : 1 : oc.num.constraints.path
    
    fprintf(fid,['\t\tcase {',int2str(ctrConstraint),'}\n']);
    
    fprintf(fid,'\t\t\tif(~isa(YL,''sym''))\n');
    fprintf(fid,['\t\t\t\t[',controlNames]);
    fprintf(fid,[',d2Hdu2] = computeControlConstraint',int2str(ctrConstraint),'(xAndLambda,const,constraint,numArcs);\n']);
    fprintf(fid,'\t\t\telse\n\t');
    fprintf(fid,['\t\t\t\t[',controlNames,',d2Hdu2] = ']);
    fprintf(fid,['deal(',controlSymNames,',0);\n']);
    fprintf(fid,'\t\t\tend\n');
    fprintf(fid,'\n');
end

fprintf(fid,'\t\tend\n');
% Output for debugging
if in.verbose
    fprintf(fid,'\t\tif numArcs > 1\n');
    
    fprintf(fid,'\t\t\tctrArc = ctrArc\n');
    fprintf(fid,'\t\t\tctrEndpoint = ctrEndpoint\n');
    
    % 	  for ctrState = 1 : 1 : oc.num.states
    % 	  	fprintf(fid,[char(oc.state.var(ctrState)),' = ',char(oc.state.var(ctrState)),'\n']);
    % 	  end
    %
    % 	  for ctrCostate = 1 : 1 : oc.num.states
    % 	  	fprintf(fid,[char(oc.costate.var(ctrCostate)),' = ',char(oc.costate.var(ctrCostate)),'\n']);
    % 	  end
    %
    for ctrControl = 1 : 1 : length(oc.control.var)
        fprintf(fid,['\t\t\t',char(oc.control.var(ctrControl)),' = ',char(oc.control.var(ctrControl)),'\n']);
    end
    %
    % names = fieldnames(in.const);
    % for ctrConst = 1 : 1 : length(names)
    % 	fprintf(fid,[char(names{ctrConst}),' = ',char(names{ctrConst}),'\n']);
    % end
    
    fprintf(fid,'\t\tend\n');
end

% Obtain H(t+) and H(t-). Can use H unconstrained since S = 0 on constraint.
% Don't need free final time condition when time is converted into a state

fprintf(fid,'\t\tif ctrEndpoint == 1\n');
if in.rootSolving == 0
    fprintf(fid,['\t\t\tH_t_plus(ctrArc) = real(',char(oc.hamiltonian.unconstrained.expression),');\n']);
else
    fprintf(fid,['\t\t\tH_t_plus(ctrArc) = (',char(oc.hamiltonian.unconstrained.expression),');\n']);
end
if in.verbose
    fprintf(fid,'\t\t\tH_t_plus(ctrArc)\n');
end
fprintf(fid,'\t\telseif ctrEndpoint == 2\n');
if in.rootSolving == 0
    fprintf(fid,['\t\t\tH_t_minus(ctrArc) = real(',char(oc.hamiltonian.unconstrained.expression),');\n']);
else
    fprintf(fid,['\t\t\tH_t_minus(ctrArc) = (',char(oc.hamiltonian.unconstrained.expression),');\n']);
end
if in.verbose
    fprintf(fid,'\t\t\tH_t_minus(ctrArc)\n');
end
fprintf(fid,'\t\tend\n');

% Information across arcs
fprintf(fid,'\t\tif ctrArc > 1 && ctrEndpoint == 1\n'); % do calculations only once for each joining arc

% Loop through all possible interior point constraints
fprintf(fid,['\t\t\tfor ctrInteriorPoint = 1 : 1 : ',int2str(oc.num.constraints.path+oc.num.constraints.interiorPoint),'\n']);

fprintf(fid,'\t\t\t\tswitch interiorPointConstraintSequence(ctrInteriorPoint,ctrArc)\n');

fprintf(fid,['\t\t\t\tcase {0}\n']);
fprintf(fid,['\t\t\t\t\tdNdx = zeros(',int2str(oc.num.states),',1);\n']);
fprintf(fid,'\t\t\t\t\tdNdt = 0;\n');
fprintf(fid,'\t\t\t\t\tpii = 0;\n');
% fprintf(fid,'\t\t\t\t\tctrArc;\n');
% for ctrConstraintDeriv1 = 1 : 1 : length(oc.constraint.interiorPoint.expression{1,1}(:,1))
%     fprintf(fid,'\t\t\t\t\tinteriorPointConstraintIndex1 = interiorPointConstraintIndex1 + 1;\n');
%     fprintf(fid,['\t\t\t\t\tconstraintDeriv1(interiorPointConstraintIndex1,1) = ',char(oc.constraint.interiorPoint.expression{1,1}(ctrConstraintDeriv1,1)),';\n']);
% end
% 
% % Compute partial of interior point constraint with respect to the state
% fprintf(fid,['\t\t\t\t\tdNdx = interiorPoint',int2str(1),'statePartial(']);
% for ctrControl = 1 : 1 : oc.num.controls
%     fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
% end
% for ctrState = 1 : 1 : oc.num.states
%     fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
%     if ctrState ~= oc.num.states
%         fprintf(fid,',');
%     end
% end
% fprintf(fid,',const,constraint);\n');
% fprintf(fid,'\n');
% % Compute partial of interior point constraint with respect to the independent variable
% fprintf(fid,['\t\t\t\t\tdNdt = interiorPoint',int2str(1),'independentVariablePartial(']);
% for ctrControl = 1 : 1 : oc.num.controls
%     fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
% end
% for ctrState = 1 : 1 : oc.num.states
%     fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
%     if ctrState ~= oc.num.states
%         fprintf(fid,',');
%     end
% end
% fprintf(fid,',const,constraint);\n');
% fprintf(fid,'\n');
% 
% % Select appropriate interior point Lagrange multipliers
% fprintf(fid,'\t\t\t\t\t\tpii1 = lmInteriorPoint1(lmInteriorPointIndex1:lmInteriorPointIndex1+1);\n');
% % fprintf(fid,'\t\t\t\t\t\tlmInteriorPointIndex1 = lmInteriorPointIndex1+interiorPointNumLagrangeMultipliers(interiorPointConstraintSequence(ctrInteriorPoint,ctrArc));\n');
% 
% %   if oc.cost.interiorPoint.constraintIndex > 0
% %     if ctrInteriorPointConstraint == (oc.num.constraints.path+oc.cost.interiorPoint.constraintIndex)
% %     	% Compute partial of interior point cost with respect to the state
% %     	fprintf(fid,['\t\t\t\t\tdIdx = interiorPointCostStatePartial(']);
% %     	for ctrControl = 1 : 1 : oc.num.controls
% %     		fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
% %     	end
% %     	for ctrState = 1 : 1 : oc.num.states
% %     		fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
% %     		if ctrState ~= oc.num.states
% %     			fprintf(fid,',');
% %     		end
% %     	end
% %     	fprintf(fid,',const,constraint);\n');
% %     	fprintf(fid,'\n');
% %     	% Compute partial of interior point cost with respect to the independent variable
% %     	fprintf(fid,['\t\t\t\t\tdIdt = interiorPointCostIndependentVariablePartial(']);
% %     	for ctrControl = 1 : 1 : oc.num.controls
% %     		fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
% %     	end
% %     	for ctrState = 1 : 1 : oc.num.states
% %     		fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
% %     		if ctrState ~= oc.num.states
% %     			fprintf(fid,',');
% %     		end
% %     	end
% %     	fprintf(fid,',const,constraint);\n');
% %     	fprintf(fid,'\n');
% %
% %       fprintf(fid,'\t\t\t\thamiltonianDiscontinuity = hamiltonianDiscontinuity + dIdt;\n');
% %       fprintf(fid,'\t\t\t\tcostateDiscontinuity = costateDiscontinuity + dIdx;\n');
% %     end
% %   end
% % end
% fprintf(fid,'\t\t\t\thamiltonianDiscontinuity = hamiltonianDiscontinuity + pii1''*dNdt;\n');
% fprintf(fid,'\t\t\t\tcostateDiscontinuity = costateDiscontinuity + dNdx*pii1;\n');
for ctrInteriorPointConstraint = 1 : 1 : oc.num.constraints.path+oc.num.constraints.interiorPoint
    fprintf(fid,['\t\t\t\tcase {',int2str(ctrInteriorPointConstraint),'}\n']);
    % Compute constraint derivative
    for ctrConstraintDeriv = 1 : 1 : length(oc.constraint.interiorPoint.expression{ctrInteriorPointConstraint,1}(:,1))
        fprintf(fid,'\t\t\t\t\tinteriorPointConstraintIndex = interiorPointConstraintIndex + 1;\n');
        fprintf(fid,['\t\t\t\t\tconstraintDeriv(interiorPointConstraintIndex,1) = ',char(oc.constraint.interiorPoint.expression{ctrInteriorPointConstraint,1}(ctrConstraintDeriv,1)),';\n']);
    end
    
    % Compute partial of interior point constraint with respect to the state
    fprintf(fid,['\t\t\t\t\tdNdx = interiorPoint',int2str(ctrInteriorPointConstraint),'statePartial(']);
    for ctrControl = 1 : 1 : oc.num.controls
        fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
    end
    for ctrState = 1 : 1 : oc.num.states
        fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
        if ctrState ~= oc.num.states
            fprintf(fid,',');
        end
    end
    fprintf(fid,',const,constraint);\n');
    fprintf(fid,'\n');
    % Compute partial of interior point constraint with respect to the independent variable
    fprintf(fid,['\t\t\t\t\tdNdt = interiorPoint',int2str(ctrInteriorPointConstraint),'independentVariablePartial(']);
    for ctrControl = 1 : 1 : oc.num.controls
        fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
    end
    for ctrState = 1 : 1 : oc.num.states
        fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
        if ctrState ~= oc.num.states
            fprintf(fid,',');
        end
    end
    fprintf(fid,',const,constraint);\n');
    fprintf(fid,'\n');
    
    % Select appropriate interior point Lagrange multipliers
    fprintf(fid,'\t\t\t\t\t\tpii = lmInteriorPoint(lmInteriorPointIndex:lmInteriorPointIndex+interiorPointNumLagrangeMultipliers(interiorPointConstraintSequence(ctrInteriorPoint,ctrArc))-1);\n');
    fprintf(fid,'\t\t\t\t\t\tlmInteriorPointIndex = lmInteriorPointIndex+interiorPointNumLagrangeMultipliers(interiorPointConstraintSequence(ctrInteriorPoint,ctrArc));\n');
    
    if oc.cost.interiorPoint.constraintIndex > 0
        if ctrInteriorPointConstraint == (oc.num.constraints.path+oc.cost.interiorPoint.constraintIndex)
            % Compute partial of interior point cost with respect to the state
            fprintf(fid,['\t\t\t\t\tdIdx = interiorPointCostStatePartial(']);
            for ctrControl = 1 : 1 : oc.num.controls
                fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
            end
            for ctrState = 1 : 1 : oc.num.states
                fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
                if ctrState ~= oc.num.states
                    fprintf(fid,',');
                end
            end
            fprintf(fid,',const,constraint);\n');
            fprintf(fid,'\n');
            % Compute partial of interior point cost with respect to the independent variable
            fprintf(fid,['\t\t\t\t\tdIdt = interiorPointCostIndependentVariablePartial(']);
            for ctrControl = 1 : 1 : oc.num.controls
                fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
            end
            for ctrState = 1 : 1 : oc.num.states
                fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
                if ctrState ~= oc.num.states
                    fprintf(fid,',');
                end
            end
            fprintf(fid,',const,constraint);\n');
            fprintf(fid,'\n');
            
            fprintf(fid,'\t\t\t\thamiltonianDiscontinuity = hamiltonianDiscontinuity + dIdt;\n');
            fprintf(fid,'\t\t\t\tcostateDiscontinuity = costateDiscontinuity + dIdx;\n');
        end
    end
end
% fprintf(fid,'\t\t\t\thamiltonianDiscontinuity = hamiltonianDiscontinuity + pii''*dNdt;\n');
% fprintf(fid,'\t\t\t\tcostateDiscontinuity = costateDiscontinuity + dNdx*pii;\n');
fprintf(fid,'\t\t\t\tend\n'); % switch interiorPointConstraintSequence

%% Change made in March 2015: size mismatch issues fixed for pii, dNdx and dNdt
fprintf(fid,'\t\t\t\thamiltonianDiscontinuity = hamiltonianDiscontinuity + pii''*dNdt;\n');
fprintf(fid,'\t\t\t\tcostateDiscontinuity = costateDiscontinuity + dNdx*pii;\n');
% fprintf(fid,'\t\t\t\thamiltonianDiscontinuity = hamiltonianDiscontinuity;\n');
% fprintf(fid,'\t\t\t\tcostateDiscontinuity = costateDiscontinuity;\n');
fprintf(fid,'\t\t\tend\n'); % ctrInteriorPoint

% Continuity of states between arcs
if ~in.convertParametersToStates
    fprintf(fid,['\t\t\tcontinuityStates(indexContinuityStates:indexContinuityStates+',int2str(oc.num.states),'-1,1) = YR(1:',int2str(oc.num.states),',ctrArc-1) - YL(1:',int2str(oc.num.states),',ctrArc);\n']);
    fprintf(fid,['\t\t\tindexContinuityStates = indexContinuityStates + ',int2str(oc.num.states),';\n']);
else
    fprintf(fid,['\t\t\tcontinuityStates(indexContinuityStates:indexContinuityStates+',int2str(oc.num.states),'+numArcs-1,1) = YR([1:',int2str(oc.num.states),',',int2str(2*oc.num.states),'+1:',int2str(2*oc.num.states),'+numArcs],ctrArc-1) - YL([1:',int2str(oc.num.states),',',int2str(2*oc.num.states),'+1:',int2str(2*oc.num.states),'+numArcs],ctrArc);\n']);
    fprintf(fid,['\t\t\tindexContinuityStates = indexContinuityStates + ',int2str(oc.num.states),'+numArcs;\n']);
end

% Discontinuity of Hamiltonian and costates across arcs
indRange = [int2str(oc.num.states+1),':',int2str(2*oc.num.states)];

fprintf(fid,'\t\t\tindexHarcs = indexHarcs + 1;\n');
fprintf(fid,'\t\t\tHarcs(indexHarcs) = H_t_minus(ctrArc-1) - H_t_plus(ctrArc) - hamiltonianDiscontinuity;\n');
fprintf(fid,['\t\t\tcostateArcs(indexCostateArcs:indexCostateArcs+',int2str(oc.num.states),'-1,1) = YR(',indRange,',ctrArc-1) - YL(',indRange,',ctrArc) - costateDiscontinuity;\n']);
fprintf(fid,['\t\t\tindexCostateArcs = indexCostateArcs + ',int2str(oc.num.states),';\n']);

fprintf(fid,'\t\tend\n'); % information across arcs

fprintf(fid,'\tend\n'); % loop ctrEndpoint

fprintf(fid,'end\n'); % loop ctrArc

%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero Vector Output %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Open vector for writing
fprintf(fid,'\n');
fprintf(fid,'if numArcs < 2\n');
if in.rootSolving == 0
    fprintf(fid,'\tzeroVec = real([');
else
    fprintf(fid,'\tzeroVec = ([');
end
% Write boundary conditions
fprintf(fid,'xAndLambda0Constraint;\n'); % Terminal state and costate constraint
fprintf(fid,'\t\txAndLambdaFConstraint;\n'); % Terminal state and costate constraint
fprintf(fid,'\t\tH_t_minus(end);\n'); % Free final time condition

if oc.num.constraints.path+oc.num.constraints.interiorPoint > 0
    fprintf(fid,'\t\tconstraintDeriv;\n'); % Constraint tangency conditions
    % fprintf(fid,'\t\tconstraintDeriv1;\n'); % Constraint tangency conditions
end
if in.maxNumArcs > 1
    fprintf(fid,'\t\tcontinuityStates;\n'); % Continuity of states between arcs
    fprintf(fid,'\t\tcostateArcs;\n'); % Continuity (or discontinuity) of costates across arcs
end
% Close vector
fprintf(fid,'\t]);\n');

fprintf(fid,'else\n');

if in.rootSolving == 0
    fprintf(fid,'\tzeroVec = real([');
else
    fprintf(fid,'\tzeroVec = ([');
end
% Write boundary conditions
fprintf(fid,'xAndLambda0Constraint;\n'); % Terminal state and costate constraint
fprintf(fid,'\t\txAndLambdaFConstraint;\n'); % Terminal state and costate constraint
fprintf(fid,'\t\tH_t_minus(end);\n'); % Free final time condition
fprintf(fid,'\t\tHarcs;\n'); % Hamiltonian continuity conditions
if oc.num.constraints.path+oc.num.constraints.interiorPoint > 0
    fprintf(fid,'\t\tconstraintDeriv;\n'); % Constraint tangency conditions
    % fprintf(fid,'\t\tconstraintDeriv1;\n'); % Constraint tangency conditions
end
if in.maxNumArcs > 1 
    fprintf(fid,'\t\tcontinuityStates;\n'); % Continuity of states between arcs
    fprintf(fid,'\t\tcostateArcs;\n'); % Continuity (or discontinuity) of costates across arcs
end
% Close vector
fprintf(fid,'\t]);\n');
% fprintf(fid,'keyboard\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%
%% Verbose Output %%
%%%%%%%%%%%%%%%%%%%%

if in.verbose
    
    % 	fprintf(fid,'if numArcs > 1\n');
    
    % 	% States
    % 	for ctrState = 1 : 1 : oc.num.states
    % 		fprintf(fid,[char(oc.state.var(ctrState)),'\n']);
    % 	end
    % 	fprintf(fid,'\n');
    %
    % 	% Costates
    % 	for ctrCostate = 1 : 1 : oc.num.states
    % 		fprintf(fid,[char(oc.costate.var(ctrCostate)),'\n']);
    % 	end
    % 	fprintf(fid,'\n');
    %
    % Controls
    for ctrControl = 1 : 1 : length(oc.control.var)
        fprintf(fid,[char(oc.control.var(ctrControl)),'\n']);
    end
    
    % 	fprintf(fid,'p\n');
    % 	% fprintf(fid,'p(end)\n');
    %
    % 	fprintf(fid,'YL\n');
    % 	fprintf(fid,'YR\n');
    
    %   fprintf(fid,'oW\n');
    %   fprintf(fid,'thrustMult/2*(thrust+1)/(g0*Isp)\n');
    
    %     fprintf(fid,'lamV\n');
    %     fprintf(fid,'v\n');
    %     fprintf(fid,'c\n');
    
    % 	% Zero vector
    fprintf(fid,'size(constraintDeriv)\n'); % Constraint tangency conditions
    fprintf(fid,'size(continuityStates)\n'); % Continuity of states between arcs
    fprintf(fid,'size(xAndLambda0Constraint)\n'); % Terminal state and costate constraint
    fprintf(fid,'size(xAndLambdaFConstraint)\n'); % Terminal state and costate constraint
    fprintf(fid,'size(H_t_minus(end))\n'); % Free final time condition
    fprintf(fid,'size(Harcs)\n'); % Continuity (or discontinuity) of Hamiltonian across arcs
    fprintf(fid,'size(costateArcs)\n'); % Continuity (or discontinuity) of costates across arcs
    %	fprintf(fid,'zeroVec\n');
    
    % fprintf(fid,'(YR(9)*(2*YR(1) - 2*YR(4)))/(abs(2*YR(4) - 2*YR(1))^2*abs(YR(9))^2)^(1/2)\n');
    
    %   fprintf(fid,'qDot\n');
    %   fprintf(fid,'v\n');
    %   fprintf(fid,'r\n');
    %   fprintf(fid,'re\n');
    %   fprintf(fid,'k\n');
    %   fprintf(fid,'rho0\n');
    %   fprintf(fid,'H\n');
    %   fprintf(fid,'rn\n');
    %   fprintf(fid,'interiorPointConstraintSequence\n');
    %   fprintf(fid,'lmInteriorPoint\n');
    
    % fprintf(fid,'if length(arcSequence) > 1\n');
    % fprintf(fid,'error(''stop'')\n');
    % fprintf(fid,'end\n');
    
    % fprintf(fid,'pause\n');
    
    % 	fprintf(fid,'keyboard\n');
    
    % 	fprintf(fid,'end\n');
    
end

%%%%%%%%%%%%%%%%%%%%
%% Close Out File %%
%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'return\n');
fprintf(fid,'\n');

return
