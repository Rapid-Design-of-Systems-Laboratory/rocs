function [] = writeDerivFunc(fid,in,oc,multipointBVP)

	% Need to change to variable arguments to support region variable

	%%%%%%%%%%%%%%%%%%
	%% Write Header %%
	%%%%%%%%%%%%%%%%%%

	% Open function call
	if multipointBVP
		fprintf(fid,'function [xDotOut] = derivFunc(t,X,region,p');
	else
		fprintf(fid,'function [xDotOut] = derivFunc(t,X,p');
	end

	% Close out function call
	if in.rootSolving == 0
		fprintf(fid,',const,constraint,arcSequence,interiorPointConstraintSequence,interiorPointNumLagrangeMultipliers,x0,xf) %%#ok<INUSD,INUSL>\n');
	else
		fprintf(fid,',const,constraint,arcSequence,varargin) %%#ok<INUSD,INUSL>\n');
	end

	% if in.rootSolving == 0
	if in.rootSolving == 1
		fprintf(fid,'if(~isa(X,''sym''))\n');
	end
		writeHeader(fid,in);
		fprintf(fid,'assert(isa(t, ''double''));\n');
		fprintf(fid,'assert(all(size(t)== [1 1]));\n');
		
		
		fprintf(fid,'assert(isa(X, ''double''));\n');
	
		fprintf(fid,'assert(isa(p, ''double''));\n');


		fprintf(fid,'assert(isa(arcSequence, ''double''));\n');
		
		if in.rootSolving == 0
			if ~in.convertParametersToStates
				fprintf(fid,['assert(all(size(X)== [',int2str(2*oc.num.states),' 1]));\n']);
			else
				fprintf(fid,['coder.varsize(''X'',[',int2str(2*oc.num.states+in.maxNumArcs),' 1],[1 0]);\n']);
			end
			fprintf(fid,['coder.varsize(''arcSequence'', [1 ',int2str(in.maxNumArcs),']);\n']);
			fprintf(fid,['coder.varsize(''p'', [',int2str(in.maxNumArcs + in.maxNumArcs*oc.num.states + oc.num.constraints.initial + oc.num.constraints.terminal),' 1]);\n']);
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
			fprintf(fid,'assert(isa(x0, ''double''));\n');
			fprintf(fid,['assert(all(size(x0)== [',int2str(oc.num.states),' 1]));\n']);
			fprintf(fid,'assert(isa(xf, ''double''));\n');
			fprintf(fid,['assert(all(size(xf)== [',int2str(oc.num.states),' 1]));\n']);
		
			if ~in.convertParametersToStates
				fprintf(fid,['xDot = nan(',num2str(2*oc.num.states),',1);\n']);
			else
				fprintf(fid,['xDot = nan(',num2str(2*oc.num.states+in.maxNumArcs),',1);\n']);
				% fprintf(fid,['assert(isa(xDot, ''double''));\n']);
				% fprintf(fid,['coder.varsize(''xDot'',[',num2str(2*oc.num.states+in.maxNumArcs),',1],[1 0]);\n']);
			end
		end

		if multipointBVP
			fprintf(fid,'assert(isa(region, ''double''));\n');
			fprintf(fid,'assert(all(size(region)== [1 1]));\n');
		else
			fprintf(fid,'region = 1;\n');
		end
		fprintf(fid,'\txDot = zeros(size(X));\n');
		if in.rootSolving == 1
			fprintf(fid,'end\n');
		end
	% else

		% fprintf(fid,['if nargin<1\n','xDot = ',int2str(2*oc.num.states),';\nreturn;\nelse\n']);
		% fprintf(fid,'if(~isa(X,''sym''))\n\txDot = zeros(size(X));\nend\n');

	% end
	fprintf(fid,['coder.varsize(''xAndLambda'',[',int2str(oc.num.states*2+in.maxNumArcs),',1],[1,0]);\n']);
	fprintf(fid,'\n');
	
	if in.verbose
		fprintf(fid,'coder.extrinsic(''keyboard'');\n');
	end

	%%%%%%%%%%%%%%%%%%%%%%
	%% Write Parameters %%
	%%%%%%%%%%%%%%%%%%%%%%

	% if ~in.convertParametersToStates

	fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
	fprintf(fid,'%%%% Parameters %%%%\n');
	fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
	fprintf(fid,'\n');

	% Time points
	fprintf(fid,'numArcs = length(arcSequence);\n');
	% for ctr = 1 : 1 : oc.numArcs

	% end
	% pIndex = oc.numArcs;
	fprintf(fid,'\n');

	% end

	fprintf(fid,'%% States\n');
	for ctrState = 1 : 1 : oc.num.states
		fprintf(fid,[char(oc.state.var(ctrState)),' = X(',int2str(ctrState),');\n']);
	end	
	fprintf(fid,'\n');

	fprintf(fid,'%% Costates\n');
	for ctrCostate = 1 : 1 : oc.num.states
		fprintf(fid,[char(oc.costate.var(ctrCostate)),' = X(',int2str(ctrCostate+oc.num.states),');\n']);
	end
	fprintf(fid,'\n');

	if in.convertParametersToStates
		fprintf(fid,'%% Independant Variables\n');
		fprintf(fid,'if(~isa(X,''sym''))\n');
		fprintf(fid,['\ttSet = NaN(',int2str(in.maxNumArcs),',1);\n']);;
		fprintf(fid,['\tfor ctr = 1 : 1 : numArcs\n']);
		fprintf(fid,['\t\ttSet(ctr) = X(',int2str(2*oc.num.states),'+ctr);\n']);
		fprintf(fid,['\tend\n']);
		fprintf(fid,['end\n']);
		fprintf(fid,'\n');
	end

	%%%%%%%%%%%%%%%%%%%%%
	%% Write Constants %%
	%%%%%%%%%%%%%%%%%%%%%
	
	if isfield(in,'const')
		writeConstants(fid,in.const,false);
		fprintf(fid,'\n');
	end
	
	if isfield(in,'constraintVal')
		writeConstraints(fid,in.constraintVal,false);
	end
	fprintf(fid,'\n');
	
	if isfield(in,'quantity')
		writeQuantities(fid,in.quantity);
		fprintf(fid,'\n');
	end
	fprintf(fid,'\n');

	if ~in.convertParametersToStates
		fprintf(fid,['tSet = p(1:numArcs);\n']);
	else
		fprintf(fid,'if(~isa(X,''sym''))\n');
		fprintf(fid,['\tt_ = abs(tSet(region));\n']);
		fprintf(fid,['else\n\tt_=sym(''t_'');\nend\n']);
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Information from Arcs %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Determine current arc index
	fprintf(fid,'indexArc = arcSequence(region);\n');

	
	% Initialize control
	for ctrControl = 1 : 1 : length(oc.control.var)
		fprintf(fid,[char(oc.control.var(ctrControl)),' = NaN;\n']);
	end

	% Assign costate derivatives based on trajectory arc
	fprintf(fid,'switch indexArc\n');

	%%%%%%%%%%%%%%%%%%%%%%%
	%% Unconstrained Arc %%
	%%%%%%%%%%%%%%%%%%%%%%%

	% Write case statement
	fprintf(fid,'\tcase {0} %% unconstrained arc\n');

	% Write control
	fprintf(fid,'\t\txAndLambda = [');
	for ctrState = 1 : 1 : oc.num.states
		fprintf(fid,[char(oc.state.var(ctrState)),';']);
	end
	for ctrCostate = 1 : 1 : oc.num.states
		fprintf(fid,char(oc.costate.var(ctrCostate)));
		if ctrCostate ~= oc.num.states
			fprintf(fid,';');
		else
			if in.convertParametersToStates
				% fprintf(fid,';NaN(length(arcSequence))');
				for ctr = 1 : 1 : in.maxNumArcs
					fprintf(fid,';NaN');
				end
			end
			fprintf(fid,'];\n');
		end
	end
	if in.convertParametersToStates
		fprintf(fid,'if(~isa(X,''sym''))\n');
		fprintf(fid,'\t\tfor ctr = 1 : 1 : length(arcSequence)\n');
		fprintf(fid,['\t\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
		fprintf(fid,'\t\tend\n');
		fprintf(fid,'end;\n');
	end
	
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
	fprintf(fid,'\t\tif(~isa(X,''sym''))\n');
	fprintf(fid,'\t\t\t[');
	fprintf(fid,controlNames);
	fprintf(fid,',hamiltonian] = computeControlUnconstrained(xAndLambda,const,constraint,length(arcSequence));\n');
	fprintf(fid,'\t\telse\n');
	fprintf(fid,['\t\t\t[',controlNames,',hamiltonian] = ']);
	fprintf(fid,['deal(',controlSymNames,',sym(''hamiltonian''));\n']);
	fprintf(fid,'\t\tend\n');
	fprintf(fid,'\n');

	% Write derivative of costates
	if in.rootSolving == 0
		fprintf(fid,['\t\txDot(',int2str(oc.num.states+1),':',int2str(2*oc.num.states),',1) = real([']);
	else
		fprintf(fid,['\t\txDot(',int2str(oc.num.states+1),':',int2str(2*oc.num.states),',1) = ([']);
	end
	for ctrCostateDeriv = 1 : 1 : oc.num.costates

		if ctrCostateDeriv ~= 1
			% Add spaces
			fprintf(fid,'\t\t\t\t');
		end

		% Write EOM
		fprintf(fid,char(oc.costate.rate.unconstrained(ctrCostateDeriv,1)));

		% Close out line
		if ctrCostateDeriv == oc.num.costates
			% Last entry
			fprintf(fid,']);\n');
		else
			% Middle entry
			fprintf(fid,'; ...\n');
		end
	end
	fprintf(fid,'\n');
	if in.convertParametersToStates
	%	fprintf(fid,['\t\txDot(',num2str(oc.num.states+oc.num.origStates),'+region) = -hamiltonian;\n']);
	end
	
	%%%%%%%%%%%%%%%%%%%%%%
	%% Constrained Arcs %%
	%%%%%%%%%%%%%%%%%%%%%%

	% Loop through each constraint type
	for ctrPathConstraint = 1 : 1 : oc.num.constraints.path

		% Write case statement for constraint type
		fprintf(fid,['\tcase {',int2str(ctrPathConstraint),'} %% constrained arc\n']);

		% Write control values
		fprintf(fid,'\t\txAndLambda = [');
		for ctrState = 1 : 1 : oc.num.states
			fprintf(fid,[char(oc.state.var(ctrState)),';']);
		end
		for ctrCostate = 1 : 1 : oc.num.states
			fprintf(fid,char(oc.costate.var(ctrCostate)));
			if ctrCostate ~= oc.num.states
				fprintf(fid,';');
			else
				if in.convertParametersToStates
					% fprintf(fid,';NaN(length(arcSequence))');
					for ctr = 1 : 1 : in.maxNumArcs
						fprintf(fid,';NaN');
					end
				end
				fprintf(fid,'];\n');
			end
		end
		if in.convertParametersToStates
			fprintf(fid,'\t\t\tif(~isa(X,''sym''))\n');
			fprintf(fid,'\t\t\t\tfor ctr = 1 : 1 : length(arcSequence)\n');
			fprintf(fid,['\t\t\t\txAndLambda(',num2str(2*oc.num.states),'+ctr) = tSet(ctr);\n']);
			fprintf(fid,'\t\t\t\tend\n');
			fprintf(fid,'\t\t\tend\n');
		end
		
		fprintf(fid,'\t\t\tif(~isa(X,''sym''))\n');
		fprintf(fid,['\t\t\t\t[',controlNames]);
		fprintf(fid,[',hamiltonian] = computeControlConstraint',int2str(ctrPathConstraint),'(xAndLambda,const,constraint,length(arcSequence));\n']);
		fprintf(fid,'\t\t\telse\n\t');
		fprintf(fid,['\t\t\t[',controlNames,',d2Hdu2] = ']);
		fprintf(fid,['deal(',controlSymNames,',0);\n']);
		fprintf(fid,'\t\t\tend\n');
		fprintf(fid,'\n');
		
		% Write muu values
		muConstraint = sym();
		fprintf(fid,['    mu',int2str(ctrPathConstraint),' = ',char(oc.lagrangeMultiplier.constraint.path{ctrPathConstraint}),';\n']);

		% Write costate derivatives
		if in.rootSolving == 0
			fprintf(fid,['\txDot(',int2str(oc.num.states+1),':',int2str(2*oc.num.states),',1) = real([']);
		else
			fprintf(fid,['\txDot(',int2str(oc.num.states+1),':',int2str(2*oc.num.states),',1) = ([']);
		end
		for ctrCostate = 1 : 1 : length(oc.costate.rate.constraint.path{ctrPathConstraint})

			if ctrCostate ~= 1
				% Add spaces
				fprintf(fid,'\t\t\t\t');
			end

			% Write EOM
			fprintf(fid,char(oc.costate.rate.constraint.path{ctrPathConstraint}(ctrCostate)));
				
			% Close out line
			if ctrCostate == length(oc.costate.rate.constraint.path{ctrPathConstraint})
				% Last entry
				fprintf(fid,']);\n');
			else
				% Middle entry
				fprintf(fid,'; ...\n');
			end
		end
		fprintf(fid,'\n');
	end

	% Close switch statement
	fprintf(fid,'end\n');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Write Equations of Motion %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	fprintf(fid,'%% Equations of motion\n');
	if in.rootSolving == 0
		fprintf(fid,['xDot(1:',int2str(oc.num.states),',1) = real([']);
	else
		fprintf(fid,['xDot(1:',int2str(oc.num.states),',1) = ([']);
	end
	for ctrStates = 1 : 1 : oc.num.states
		if ctrStates ~= 1
			% Add spaces
			fprintf(fid,'\t\t\t\t');
		end

		% Write EOM
		fprintf(fid,char(oc.state.rate(ctrStates)));

		% Close out line
		if ctrStates == oc.num.states
			fprintf(fid,']);\n');
			% Last entry
			if in.convertParametersToStates
				fprintf(fid,'if(~isa(X,''sym''))\n\tfor ctr = 1 : 1 : length(arcSequence)\n');
				fprintf(fid,['\t\txDot(',num2str(2*oc.num.states),'+ctr) = 0;\n']);
				fprintf(fid,'\tend\nend\n');
			end			
		else
			% Middle entry
			fprintf(fid,'; ...\n');
		end

	end
	fprintf(fid,'\n');
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Write "equations of motion" for the time parameters %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Multiply Derivatives By Segment Time Since Performing Fixed Time Integration %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	fprintf(fid,'%% Account for variable endpoints in derivative\n');
	if ~in.convertParametersToStates
		fprintf(fid,'xDotOut = xDot*tSet(region);\n');
	else
		% fprintf(fid,['xDot(1:',num2str(oc.num.states+oc.num.origStates),') = xDot(1:',num2str(oc.num.states+oc.num.origStates),')*t_;\n']);
		fprintf(fid,'if(~isa(X,''sym''))\n');
		fprintf(fid,['\txDotOut = xDot(1:',num2str(2*oc.num.states),'+length(arcSequence))*t_;\n']);
		fprintf(fid,'else\n');
		fprintf(fid,['\txDotOut = xDot(1:',num2str(2*oc.num.states),')*t_;\n']);
		fprintf(fid,'end;\n');
	end
	fprintf(fid,'\n');
	
	%%%%%%%%%%%%%%%%%%%%
	%% Verbose Output %%
	%%%%%%%%%%%%%%%%%%%%
	
	if in.verbose
		
		% fprintf(fid,'if max(abs(xDotOut)) > 100\n');
		% fprintf(fid,'  X\n');
		% fprintf(fid,'  xDotOut\n');
		% fprintf(fid,'  keyboard\n');
		% fprintf(fid,'end\n');
		
  end
  
%   fprintf(fid,'if gam < -88*pi/180\n');
%   fprintf(fid,'  gam*180/pi\n');
%   fprintf(fid,'end\n');

	%%%%%%%%%%%%%%%%%%%%
	%% Close Out File %%
	%%%%%%%%%%%%%%%%%%%%

	fprintf(fid,'return\n');
	if in.rootSolving == 1
		fprintf(fid,'end\n');
	end
	fprintf(fid,'\n');

	return

