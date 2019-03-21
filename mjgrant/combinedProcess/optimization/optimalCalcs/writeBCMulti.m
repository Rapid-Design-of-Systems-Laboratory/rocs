function [] = writeBC(fid,in,oc)

%%%%%%%%%%%%%%%%%%
%% Write Header %%
%%%%%%%%%%%%%%%%%%

if in.convertParametersToStates
	fprintf(fid,'function [zeroVec] = bcMulti(YL,YR,p,const,constraint,arcSequence,interiorPointConstraintSequence,interiorPointNumLagrangeMultipliers,x0,xf,ctrArc)\n');
else
	fprintf(fid,'function [zeroVec] = bcMulti(YL,YR,p,const,constraint,arcSequence,interiorPointConstraintSequence,interiorPointNumLagrangeMultipliers,x0,xf,ctrArc)\n');
end

if in.rootSolving == 0
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
	fprintf(fid,'assert(isa(interiorPointNumLagrangeMultipliers, ''double''));\n');
	if oc.num.constraints.path+oc.num.constraints.interiorPoint > 0
		fprintf(fid,['coder.varsize(''interiorPointConstraintSequence'', [',int2str(oc.num.constraints.path+oc.num.constraints.interiorPoint),' ',int2str(in.maxNumArcs),']);\n']);
	else
		fprintf(fid,['assert(all(size(interiorPointConstraintSequence) == [0 0]));\n']);
	end
	if length(oc.num.lagrangeMultipliers.interiorPoint) > 0
		fprintf(fid,['coder.varsize(''interiorPointNumLagrangeMultipliers'', [1 ',int2str(length(oc.num.lagrangeMultipliers.interiorPoint)),']);\n']);
	else
		fprintf(fid,['assert(all(size(interiorPointNumLagrangeMultipliers) == [0 0]));\n']);
	end
	
	
	fprintf(fid,['coder.varsize(''constraintDeriv'', [',int2str(in.maxNumArcs*oc.num.states),' 1]);\n']);
	% if ~in.convertParametersToStates
		fprintf(fid,['coder.varsize(''H_t_plus'', [',int2str(in.maxNumArcs),' 1]);\n']);
		fprintf(fid,['coder.varsize(''H_t_minus'', [',int2str(in.maxNumArcs),' 1]);\n']);
		fprintf(fid,['coder.varsize(''Harcs'', [',int2str(in.maxNumArcs),' 1]);\n']);
	% end

	fprintf(fid,['coder.varsize(''continuityStates'', [',int2str(in.maxNumArcs*oc.num.states),' 1]);\n']);
	fprintf(fid,['coder.varsize(''costateArcs'', [',int2str(in.maxNumArcs*oc.num.states),' 1]);\n']);
	fprintf(fid,['coder.varsize(''dNdt'', [',int2str(in.maxNumArcs*oc.num.states),',1]);\n']);
	fprintf(fid,['coder.varsize(''dNdx'', [',int2str(in.maxNumArcs*oc.num.states),',',int2str(oc.num.states),']);\n']);
	
	fprintf(fid,['coder.varsize(''hamiltonianDiscontinuity'', [',int2str(in.maxNumArcs*oc.num.states),',1]);\n']);
	fprintf(fid,['coder.varsize(''costateDiscontinuity'', [',int2str(in.maxNumArcs*oc.num.states),',',int2str(oc.num.states),']);\n']);
	
	fprintf(fid,['assert(isa(ctrArc,''double''));\n']);
	if in.verbose
		fprintf(fid,'coder.extrinsic(''keyboard'');\n');
	end
	
end

fprintf(fid,'\n');

%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Number of arcs
fprintf(fid,'numArcs = length(arcSequence);\n');

fprintf(fid,'interiorPointConstraintIndex = 0;\n');
fprintf(fid,'lmInteriorPointIndex = 1;\n');
	
fprintf(fid,'indexConstraintDeriv = 0;\n'); % used to keep track of number of tangency constraints
fprintf(fid,'indexContinuityStates = 1;\n'); % used to keep track of number of state continuity equations across arcs
fprintf(fid,'indexCostateArcs = 1;\n'); % used to keep track of costate conditions across arcs

fprintf(fid,'indexHarcs = 0;\n'); % used to keep track of Hamiltonian conditions across arcs
fprintf(fid,'H_t_plus = NaN(numArcs,1);\n');
fprintf(fid,'H_t_minus = NaN(numArcs,1);\n');

if ~in.convertParametersToStates
	fprintf(fid,['continuityStates = NaN((numArcs-1)*',int2str(oc.num.states),',1);\n']);
else
	fprintf(fid,['continuityStates = NaN((numArcs-1)*(',int2str(oc.num.states),'+numArcs),1);\n']);
end
fprintf(fid,['costateArcs = NaN((numArcs-1)*',int2str(oc.num.states),',1);\n']);	

if in.maxNumArcs > 1
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

fprintf(fid,['xAndLambda0Constraint = NaN(',int2str(oc.num.constraints.initial+oc.num.states),',1);\n']);
fprintf(fid,['xAndLambdaFConstraint = NaN(',int2str(oc.num.constraints.terminal+oc.num.states),',1);\n']);

if ~in.convertParametersToStates
	fprintf(fid,['xAndLambda = NaN(',int2str(2*oc.num.states),',1);\n']);
else
	fprintf(fid,['xAndLambda = NaN(',int2str(2*oc.num.states+in.maxNumArcs),',1);\n']);
end
fprintf(fid,'pii = NaN;\n');

%%%%%%%%%%%%%%%%%%%%%%
%% Write Parameters %%
%%%%%%%%%%%%%%%%%%%%%%
	
	fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
	fprintf(fid,'%%%% Parameters %%%%\n');
	fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
	fprintf(fid,'\n');
	
	% Time points
	% if ~in.convertParametersToStates
	% 	fprintf(fid,['tSet = p(1:numArcs);\n']);
	% else
	% 	fprintf(fid,['t_ = X(',num2str(oc.num.origStates),'+region,1);']);
	% end
	% fprintf(fid,'\n');
	
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
		fprintf(fid,['lmInteriorPoint = p(numArcs+',int2str(oc.num.constraints.initial),'+',int2str(oc.num.constraints.terminal),'+1:end);\n']);
	else
		fprintf(fid,['lmInteriorPoint = p(0+',int2str(oc.num.constraints.initial),'+',int2str(oc.num.constraints.terminal),'+1:end);\n']);
	end
	fprintf(fid,'\n');
	
	fprintf(fid,'dNdt = NaN(length(lmInteriorPoint),1);\n');
	fprintf(fid,['dNdx = NaN(',int2str(oc.num.states),',length(lmInteriorPoint));\n']);
	fprintf(fid,'constraintDeriv = NaN(length(lmInteriorPoint),1);\n');
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write Constants and Constraint Values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writeConstants(fid,in.const,false);
fprintf(fid,'\n');
fprintf(fid,'\n');
if isfield(in,'constraintVal')
	writeConstraints(fid,in.constraintVal,false);
end
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Through Each Arc and Obtain Boundary Condition Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Also loop through left and right endpoint to reduce repeated code
fprintf(fid,'if ctrArc == 1\n');
	
fprintf(fid,'\t%% Left endpoint\n');
fprintf(fid,'\t%% States\n');
for ctrState = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.state.var(ctrState)),' = YL(',int2str(ctrState),',ctrArc);\n']);
end
fprintf(fid,'\n');

fprintf(fid,'\t%% Costates\n');
for ctrCostate = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.costate.var(ctrCostate)),' = YL(',int2str(ctrCostate+oc.num.states),',ctrArc);\n']);
end
fprintf(fid,'\n');

if in.convertParametersToStates
	fprintf(fid,'\t%% Independant Variables\n');
	fprintf(fid,'\tif(~isa(YL,''sym''))\n');
	fprintf(fid,['\t\ttSet = NaN(',int2str(in.maxNumArcs),',1);\n']);
	fprintf(fid,'\tend;\n');
	fprintf(fid,['\tfor ctr = 1 : 1 : numArcs\n']);
	fprintf(fid,['\t\ttSet(ctr) = YL(',int2str(2*oc.num.states),'+ctr,ctrArc);\n']);
	fprintf(fid,['\tend\n']);
	fprintf(fid,'\n');
end

% Write control input for boundary point
fprintf(fid,'\txAndLambda = [');
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
	fprintf(fid,'\tfor ctr = 1 : 1 : numArcs\n');
	fprintf(fid,['\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
	fprintf(fid,'\tend\n');
end

% Initial state and costate constraint
fprintf(fid,'\txAndLambda0Constraint = [\n');
for ctrState = 1 : 1 : oc.num.states
	if ismember(oc.state.var(ctrState),oc.constraint.initial) % Initial state specified
		fprintf(fid,['\t\t',char(oc.state.var(ctrState)),' - x0(',int2str(ctrState),');\n']);
	end
	fprintf(fid,['\t\t',char(oc.costate.var(ctrState)),' - ',char(oc.costate.initial(ctrState)),';\n']);
end

fprintf(fid,'\t];\n');
fprintf(fid,'\tctrArc = numArcs;');
fprintf(fid,'\t%% Right endpoint\n');
fprintf(fid,'\t%% States\n');
for ctrStates = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.state.var(ctrStates)),' = YR(',int2str(ctrStates),',ctrArc);\n']);
end
if in.convertParametersToStates
	for arcNum = 1:in.maxNumArcs
		fprintf(fid,['\t',char(oc.independantVariable),num2str(arcNum),' = YR(',int2str(oc.num.states+arcNum),',ctrArc);\n']);
	end
end

fprintf(fid,'\n');

fprintf(fid,'\t%% Costates\n');
for ctrCostates = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.costate.var(ctrCostates)),' = YR(',int2str(ctrCostates+oc.num.states),',ctrArc);\n']);
end
fprintf(fid,'\n');

if in.convertParametersToStates
	fprintf(fid,'\t%% Independant Variables\n');
	fprintf(fid,'\tif(~isa(YR,''sym''))\n');
	fprintf(fid,['\t\ttSet = NaN(',int2str(in.maxNumArcs),',1);\n']);
	fprintf(fid,'\tend;\n');
	fprintf(fid,['\tfor ctr = 1 : 1 : numArcs\n']);
	fprintf(fid,['\t\ttSet(ctr) = YL(',int2str(2*oc.num.states),'+ctr,ctrArc);\n']);
	fprintf(fid,['\tend\n']);
	fprintf(fid,'\n');
end

% Write control input for boundary point
fprintf(fid,'\txAndLambda = [');
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
	fprintf(fid,'\tfor ctr = 1 : 1 : numArcs\n');
	fprintf(fid,['\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
	fprintf(fid,'\tend\n');
end
% Compute control
fprintf(fid,'\t\tswitch arcSequence(ctrArc)\n');
	
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
fprintf(fid,'\t\tif(~isa(YL,''sym''))\n');
fprintf(fid,'\t\t\t[');
fprintf(fid,controlNames);
fprintf(fid,',d2Hdu2] = computeControlUnconstrained(xAndLambda,const,constraint,numArcs);\n');
fprintf(fid,'\t\telse\n\t');
fprintf(fid,['\t\t[',controlNames,',d2Hdu2] = ']);
fprintf(fid,['deal(',controlSymNames,',0);\n']);
fprintf(fid,'\t\tend\n');
fprintf(fid,'\n');

for ctrConstraint = 1 : 1 : oc.num.constraints.path

	fprintf(fid,['\tcase {',int2str(ctrConstraint),'}\n']);

	fprintf(fid,'\t\tif(~isa(YL,''sym''))\n');
	fprintf(fid,['\t\t\t[',controlNames]);
	fprintf(fid,[',d2Hdu2] = computeControlConstraint',int2str(ctrConstraint),'(xAndLambda,const,constraint,numArcs);\n']);
	fprintf(fid,'\t\telse\n\t');
	fprintf(fid,['\t\t\t[',controlNames,',d2Hdu2] = ']);
	fprintf(fid,['deal(',controlSymNames,',0);\n']);
	fprintf(fid,'\t\tend\n');
	fprintf(fid,'\n');
end

fprintf(fid,'\t\tend\n');
if in.rootSolving == 0
	fprintf(fid,['\t\t\tH_t_minus(ctrArc) = real(',char(oc.hamiltonian.unconstrained.expression),');\n']);
else
	fprintf(fid,['\t\t\tH_t_minus(ctrArc) = (',char(oc.hamiltonian.unconstrained.expression),');\n']);
end

% Final state and costate constraint
fprintf(fid,'\txAndLambdaFConstraint = [\n');
for ctrState = 1 : 1 : oc.num.states
	if ismember(oc.state.var(ctrState),oc.constraint.terminal) % Final state specified
		fprintf(fid,['\t\t',char(oc.state.var(ctrState)),' - xf(',int2str(ctrState),');\n']);
	end
	fprintf(fid,['\t\t',char(oc.costate.var(ctrState)),' - ',char(oc.costate.terminal(ctrState)),';\n']);
end
fprintf(fid,'\t];\n');
fprintf(fid,'\n');
fprintf(fid,'\tzeroVec = [xAndLambda0Constraint; xAndLambdaFConstraint; H_t_minus(numArcs)];\n');
fprintf(fid,'\treturn;\n');
fprintf(fid,'else\n');

fprintf(fid,'\tindexArc = arcSequence(ctrArc);\n');
fprintf(fid,'\thamiltonianDiscontinuity = zeros(1,1);\n');
fprintf(fid,['\tcostateDiscontinuity = zeros(',int2str(oc.num.states),',1);\n']);

% Obtain H(t+) and H(t-). Can use H unconstrained since S = 0 on constraint.	
fprintf(fid,'\t%% Right endpoint of previous arc for H(t-)\n');
fprintf(fid,'\t%% States\n');
for ctrStates = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.state.var(ctrStates)),' = YR(',int2str(ctrStates),',ctrArc-1);\n']);
end
if in.convertParametersToStates
	for arcNum = 1:in.maxNumArcs
		fprintf(fid,['\t',char(oc.independantVariable),num2str(arcNum),' = YR(',int2str(oc.num.states+arcNum),',ctrArc-1);\n']);
	end
end

fprintf(fid,'\n');

fprintf(fid,'\t%% Costates\n');
for ctrCostates = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.costate.var(ctrCostates)),' = YR(',int2str(ctrCostates+oc.num.states),',ctrArc-1);\n']);
end
fprintf(fid,'\n');

if in.convertParametersToStates
	fprintf(fid,'\t%% Independant Variables\n');
	fprintf(fid,'\tif(~isa(YR,''sym''))\n');
	fprintf(fid,['\t\ttSet = NaN(',int2str(in.maxNumArcs),',1);\n']);
	fprintf(fid,'\tend\n');
	fprintf(fid,['\tfor ctr = 1 : 1 : numArcs\n']);
	fprintf(fid,['\t\ttSet(ctr) = YL(',int2str(2*oc.num.states),'+ctr,ctrArc-1);\n']);
	fprintf(fid,['\tend\n']);
	fprintf(fid,'\n');
end

% Write control input for boundary point
fprintf(fid,'\txAndLambda = [');
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
	fprintf(fid,'\tfor ctr = 1 : 1 : numArcs\n');
	fprintf(fid,['\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
	fprintf(fid,'\tend\n');
end
% Compute control
fprintf(fid,'\tswitch arcSequence(ctrArc-1)\n');

fprintf(fid,'\tcase {0} %% unconstrained arc\n');

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
fprintf(fid,'\t\tif(~isa(YL,''sym''))\n');
fprintf(fid,'\t\t\t[');
fprintf(fid,controlNames);
fprintf(fid,',d2Hdu2] = computeControlUnconstrained(xAndLambda,const,constraint,numArcs);\n');
fprintf(fid,'\t\telse\n\t');
fprintf(fid,['\t\t[',controlNames,',d2Hdu2] = ']);
fprintf(fid,['deal(',controlSymNames,',0);\n']);
fprintf(fid,'\t\tend\n');
fprintf(fid,'\n');

for ctrConstraint = 1 : 1 : oc.num.constraints.path

	fprintf(fid,['\tcase {',int2str(ctrConstraint),'}\n']);

	fprintf(fid,'\tif(~isa(YL,''sym''))\n');
	fprintf(fid,['\t\t[',controlNames]);
	fprintf(fid,[',d2Hdu2] = computeControlConstraint',int2str(ctrConstraint),'(xAndLambda,const,constraint,numArcs);\n']);
	fprintf(fid,'\t\t\telse\n\t');
	fprintf(fid,['\t\t[',controlNames,',d2Hdu2] = ']);
	fprintf(fid,['deal(',controlSymNames,',0);\n']);
	fprintf(fid,'\t\tend\n');
	fprintf(fid,'\n');
end

fprintf(fid,'\tend\n');
if in.rootSolving == 0
	fprintf(fid,['\tH_t_minus(ctrArc-1) = real(',char(oc.hamiltonian.unconstrained.expression),');\n']);
else
	fprintf(fid,['\tH_t_minus(ctrArc-1) = (',char(oc.hamiltonian.unconstrained.expression),');\n']);
end
fprintf(fid,'\t% Left endpoint of current arc for H(t+)');
fprintf(fid,'\t%% States\n');
for ctrState = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.state.var(ctrState)),' = YL(',int2str(ctrState),',ctrArc);\n']);
end
fprintf(fid,'\n');

fprintf(fid,'\t%% Costates\n');
for ctrCostate = 1 : 1 : oc.num.states
	fprintf(fid,['\t',char(oc.costate.var(ctrCostate)),' = YL(',int2str(ctrCostate+oc.num.states),',ctrArc);\n']);
end
fprintf(fid,'\n');

if in.convertParametersToStates
	fprintf(fid,'\t%% Independant Variables\n');
	fprintf(fid,'\tif(~isa(YL,''sym''))\n');
	fprintf(fid,['\t\ttSet = NaN(',int2str(in.maxNumArcs),',1);\n']);
	fprintf(fid,'\tend;\n');
	fprintf(fid,['\tfor ctr = 1 : 1 : numArcs\n']);
	fprintf(fid,['\t\ttSet(ctr) = YL(',int2str(2*oc.num.states),'+ctr,ctrArc);\n']);
	fprintf(fid,['\tend\n']);
	fprintf(fid,'\n');
end

% Write control input for boundary point
fprintf(fid,'\txAndLambda = [');
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
	fprintf(fid,'\tfor ctr = 1 : 1 : numArcs\n');
	fprintf(fid,['\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
	fprintf(fid,'\tend\n');
end

% Compute control
fprintf(fid,'\tswitch indexArc\n');
	
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
fprintf(fid,'\tif(~isa(YL,''sym''))\n');
fprintf(fid,'\t\t[');
fprintf(fid,controlNames);
fprintf(fid,',d2Hdu2] = computeControlUnconstrained(xAndLambda,const,constraint,numArcs);\n');
fprintf(fid,'\telse\n\t');
fprintf(fid,['\t[',controlNames,',d2Hdu2] = ']);
fprintf(fid,['deal(',controlSymNames,',0);\n']);
fprintf(fid,'\tend\n');
fprintf(fid,'\n');

for ctrConstraint = 1 : 1 : oc.num.constraints.path

	fprintf(fid,['\t\tcase {',int2str(ctrConstraint),'}\n']);

	fprintf(fid,'\tif(~isa(YL,''sym''))\n');
	fprintf(fid,['\t\t[',controlNames]);
	fprintf(fid,[',d2Hdu2] = computeControlConstraint',int2str(ctrConstraint),'(xAndLambda,const,constraint,numArcs);\n']);
	fprintf(fid,'\telse\n\t');
	fprintf(fid,['\t\t[',controlNames,',d2Hdu2] = ']);
	fprintf(fid,['deal(',controlSymNames,',0);\n']);
	fprintf(fid,'\tend\n');
	fprintf(fid,'\n');
end
fprintf(fid,'\tend\n');

if in.rootSolving == 0
	fprintf(fid,['\tH_t_plus(ctrArc) = real(',char(oc.hamiltonian.unconstrained.expression),');\n']);
else
	fprintf(fid,['\tH_t_plus(ctrArc) = (',char(oc.hamiltonian.unconstrained.expression),');\n']);
end
	
if in.verbose
	fprintf(fid,'\tH_t_plus(ctrArc)\n');
end	
% Loop through all possible interior point constraints
fprintf(fid,['\tfor ctrInteriorPoint = 1 : 1 : ',int2str(oc.num.constraints.path+oc.num.constraints.interiorPoint),'\n']);
	
fprintf(fid,'\t\tswitch interiorPointConstraintSequence(ctrInteriorPoint,ctrArc)\n');
	
fprintf(fid,['\t\tcase {0}\n']);
fprintf(fid,['\t\t\tdNdx = zeros(',int2str(oc.num.states),',1);\n']);
fprintf(fid,'\t\t\tdNdt = 0;\n');
fprintf(fid,'\t\t\tpii = 0;\n');
	
for ctrInteriorPointConstraint = 1 : 1 : oc.num.constraints.path+oc.num.constraints.interiorPoint
	fprintf(fid,['\t\tcase {',int2str(ctrInteriorPointConstraint),'}\n']);
	% Compute constraint derivative
	for ctrConstraintDeriv = 1 : 1 : length(oc.constraint.interiorPoint.expression{ctrInteriorPointConstraint,1}(:,1))
		fprintf(fid,'\t\t\tinteriorPointConstraintIndex = interiorPointConstraintIndex + 1;\n');
		fprintf(fid,['\t\t\tconstraintDeriv(interiorPointConstraintIndex,1) = ',char(oc.constraint.interiorPoint.expression{ctrInteriorPointConstraint,1}(ctrConstraintDeriv,1)),';\n']);
	end
	% Compute partial of interior point constraint with respect to the state
	fprintf(fid,['\t\t\tdNdx = interiorPoint',int2str(ctrInteriorPointConstraint),'statePartial(']);
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
	fprintf(fid,['\t\t\tdNdt = interiorPoint',int2str(ctrInteriorPointConstraint),'independentVariablePartial(']);
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
	fprintf(fid,'\t\t\t\tpii = lmInteriorPoint(lmInteriorPointIndex:lmInteriorPointIndex+interiorPointNumLagrangeMultipliers(interiorPointConstraintSequence(ctrInteriorPoint,ctrArc))-1);\n');
	fprintf(fid,'\t\t\t\tlmInteriorPointIndex = lmInteriorPointIndex+interiorPointNumLagrangeMultipliers(interiorPointConstraintSequence(ctrInteriorPoint,ctrArc));\n');
end
fprintf(fid,'\t\tend\n'); % switch interiorPointConstraintSequence
fprintf(fid,'\t\thamiltonianDiscontinuity = hamiltonianDiscontinuity + pii*dNdt;\n');
fprintf(fid,'\t\tcostateDiscontinuity = costateDiscontinuity + dNdx*pii.'';\n');
	
fprintf(fid,'\t\t\tend\n'); % ctrInteriorPoint
	
% Continuity of states between arcs
if ~in.convertParametersToStates
	fprintf(fid,['\tcontinuityStates(indexContinuityStates:indexContinuityStates+',int2str(oc.num.states),'-1,1) = YR(1:',int2str(oc.num.states),',ctrArc-1) - YL(1:',int2str(oc.num.states),',ctrArc);\n']);
	fprintf(fid,['\tindexContinuityStates = indexContinuityStates + ',int2str(oc.num.states),';\n']);
else
	fprintf(fid,['\tcontinuityStates(indexContinuityStates:indexContinuityStates+',int2str(oc.num.states),'+numArcs-1,1) = YR([1:',int2str(oc.num.states),',',int2str(2*oc.num.states),'+1:',int2str(2*oc.num.states),'+numArcs],ctrArc-1) - YL([1:',int2str(oc.num.states),',',int2str(2*oc.num.states),'+1:',int2str(2*oc.num.states),'+numArcs],ctrArc);\n']);
	fprintf(fid,['\tindexContinuityStates = indexContinuityStates + ',int2str(oc.num.states),'+numArcs;\n']);	
end

% Discontinuity of Hamiltonian and costates across arcs
indRange = [int2str(oc.num.states+1),':',int2str(2*oc.num.states)];

fprintf(fid,'\tindexHarcs = indexHarcs + 1;\n');
fprintf(fid,'\tHarcs(indexHarcs) = H_t_minus(ctrArc-1) - H_t_plus(ctrArc) - hamiltonianDiscontinuity;\n');

fprintf(fid,['\tcostateArcs(indexCostateArcs:indexCostateArcs+',int2str(oc.num.states),'-1,1) = YR(',indRange,',ctrArc-1) - YL(',indRange,',ctrArc) - costateDiscontinuity;\n']);
fprintf(fid,['\tindexCostateArcs = indexCostateArcs + ',int2str(oc.num.states),';\n']);

%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero Vector Output %%
%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'\tzeroVec = real([');
% Write boundary conditions
fprintf(fid,'\t\tHarcs;\n'); % Hamiltonian continuity conditions
if oc.num.constraints.path+oc.num.constraints.interiorPoint > 0
	fprintf(fid,'\t\tconstraintDeriv;\n'); % Constraint tangency conditions
end
if in.maxNumArcs > 1
	fprintf(fid,'\t\tcontinuityStates;\n'); % Continuity of states between arcs
	fprintf(fid,'\t\tcostateArcs\n'); % Continuity (or discontinuity) of costates across arcs
end
% Close vector
fprintf(fid,'\t]);\n');
fprintf(fid,'\treturn\n;');
fprintf(fid,'end\n'); % branch ctrArc

%%%%%%%%%%%%%%%%%%%%
%% Verbose Output %%
%%%%%%%%%%%%%%%%%%%%

if in.verbose

	fprintf(fid,'if numArcs > 1\n');

	% States
	for ctrState = 1 : 1 : oc.num.states
		fprintf(fid,[char(oc.state.var(ctrState)),'\n']);
	end
	fprintf(fid,'\n');

	% Costates
	for ctrCostate = 1 : 1 : oc.num.states
		fprintf(fid,[char(oc.costate.var(ctrCostate)),'\n']);
	end
	fprintf(fid,'\n');
	
	% Controls
	for ctrControl = 1 : 1 : length(oc.control.var)
		fprintf(fid,[char(oc.control.var(ctrControl)),'\n']);
	end
	
	fprintf(fid,'p\n');
	fprintf(fid,'p(end)\n');
	
	% Zero vector
	fprintf(fid,'constraintDeriv\n'); % Constraint tangency conditions
	fprintf(fid,'continuityStates\n'); % Continuity of states between arcs
	fprintf(fid,'xAndLambda0Constraint\n'); % Terminal state and costate constraint
	fprintf(fid,'xAndLambdaFConstraint\n'); % Terminal state and costate constraint
	fprintf(fid,'H_t_minus(end)\n'); % Free final time condition
	fprintf(fid,'Harcs\n'); % Continuity (or discontinuity) of Hamiltonian across arcs
	fprintf(fid,'costateArcs\n'); % Continuity (or discontinuity) of costates across arcs
	fprintf(fid,'zeroVec\n');
	
	% fprintf(fid,'if length(arcSequence) > 1\n');
	% fprintf(fid,'error(''stop'')\n');
	% fprintf(fid,'end\n');
	
	% fprintf(fid,'pause\n');
	
	% fprintf(fid,'keyboard\n');
	
	fprintf(fid,'end\n');
	
end

%%%%%%%%%%%%%%%%%%%%
%% Close Out File %%
%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'return\n');
fprintf(fid,'\n');

return