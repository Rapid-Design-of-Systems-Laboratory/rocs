function [ out_str ] = makeSplitStmKernel(eom_func, in, out, varargin)
% Generates C code for the EOM file from the symbolic expressions 
%
%   Takes state equation expressions and also generates expressions for 
%   n-th order Jacobians (depending on config setting)
%
% eom_func   : Valid function handle to system dynamics file
% config     : Structure containing configuration information
%
% jac_str    : Contains the final C code output
%
syms x t;
jac_str = '';
states_str = '';

% Copy in required variables from configuration
nStatesAndCostates  = 2*in.oc.num.states;	% States + Co-states
maxSttOrder = in.gpuSolve.sttOrder;

% Create symbolic variables for the states
x = sym('x',[nStatesAndCostates 1]);

% Create symbolic variables for the auxiliary data
const_fields = fieldnames(in.const);
const_sym = struct();
for i=1:numel(const_fields)
	fieldname = char(const_fields(i));
	const_sym.(fieldname) = sym(fieldname);
end

% Obtain fieldnames of constraints
if isfield(in,'constraintVal')
	constraint_names = fieldnames(in.constraintVal);
else
	constraint_names = {};
end
constraint_sym = struct();
for i=1:numel(constraint_names)
	fieldname = char(constraint_names(i));
	constraint_sym.(fieldname) = sym(fieldname);
end

in.const_sym = const_sym;
in.const_fields = const_fields;
in.constraint_sym = constraint_sym;
in.constraint_names = constraint_names;

timeParams = sym('t',[1 in.maxNumArcs]);

% Write the constants
const_str = '';
for i=1:numel(const_fields)
	fieldname = char(const_fields(i));
	c_fieldname = sym2cexp(const_sym.(fieldname),true);
	const_str = sprintf('%s\tconst double %s = d_const[%d];\n',const_str,c_fieldname,i-1);
end
for i=1:numel(constraint_names)
	fieldname = char(constraint_names(i));
	c_fieldname = sym2cexp(constraint_sym.(fieldname),true);
	const_str = sprintf('%s\tconst double %s = d_constraints[%d];\n',const_str,c_fieldname,i-1);
end

init_str = '';
for i=1:nStatesAndCostates
	init_str = sprintf('%s\t_num_t x%d = d_data[%d];\n',init_str,i,i-1);
end

states_str = [init_str,const_str];
states_str = sprintf('%s\tconst _num_t t_ = fabs(d_data[%d+_arcId]);\n',states_str,2*out.oc.num.states);

states_str = sprintf('%s\t\t_num_t ',states_str);
for ctrControl = 1 : 1 : out.oc.num.controls
	states_str = sprintf('%s%s, ',states_str,char(out.oc.control.var(ctrControl)));
end
states_str = strcat(states_str,'hamiltonian;\n');

states_str = sprintf('%sswitch(_arcType){\n',states_str);

if ~isfield(in.oc.constraint,'path')
	arcTypes = 0;
else
	arcTypes = 0:size(in.oc.constraint.path,1);
end
for arcIndex = 1:length(arcTypes)
	arcType = arcTypes(arcIndex);
	tVars = sym(char(out.oc.independentVariable),[arcIndex 1]);
	fprintf('Arc Type %d\n',arcType);
	states_str = sprintf('%s\t\tcase %d: // Arc Type %d\n',states_str,arcType,arcType);
	
	% Compute control for the arc
	
	if arcType == 0
		controlfn = 'computeControlUnconstrained';
	else
		controlfn = ['computeControlConstraint',int2str(arcType)];
	end

	states_str = sprintf('%s\t\t\t%s(d_data,d_const,d_constraints,',states_str,controlfn);
	for ctrControl = 1 : 1 : out.oc.num.controls
		states_str = strcat(states_str,['&',char(out.oc.control.var(ctrControl)),',']);
	end
	states_str = strcat(states_str,'&hamiltonian);\n');
	
	% Evaluate the state equations
	fx = eom_func(t,[x;tVars],arcIndex,timeParams,const_sym,constraint_sym,arcTypes);
	states_str = sprintf('%s\t\tswitch(_i){ \n',states_str);
	for i=1:numel(fx)
		cexp = sym2cexp(fx(i));
		states_str = sprintf('%s\t\t\tcase %d:\n\t\t\treturn %s\n',states_str,i-1,cexp);
	end
	states_str = sprintf('%s\t\t\tdefault:\n\t\t\treturn 0.0;\n\t\t}\n',states_str);
end
states_str = strcat(states_str,'\t}\treturn 0;\n');

fprintf('Writing file to CU\n');

if strcmpi(in.gpuSolve.derivativeMethod,'csd')
	out_str = fileread('splitstm_kernel.csd.tmpl.cu');
	out_str = strrep(out_str,'%DATATYPE%','complex_t');
else
	out_str = fileread('splitstm_kernel.fd.tmpl.cu');
	out_str = strrep(out_str,'%DATATYPE%','double');
end

out_str = strrep(out_str,'%STATE_FUNCTIONS%',states_str);
out_str = strrep(out_str,'%NUMSTATES%',int2str(nStatesAndCostates));
out_str = strrep(out_str,'%MAXARCS%',int2str(in.maxNumArcs));
	

end