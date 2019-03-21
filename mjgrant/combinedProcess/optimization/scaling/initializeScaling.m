function [OUT] = initializeScaling(in,out)

% scaleVal = in.scale;
for ctrScale = 1 : 1 : length(in.scale(:,1))
	if isnumeric(in.scale{ctrScale,2}) % constant scaling factor
		% in.vars.scaleValFn{ctr1,2} = in.scale{ctr1,2};
	elseif ~isempty(strfind(in.scale{ctrScale,2},'.')) % value from structure in input file
	else % comes from trajectory states
		% Create symbolic function
		in.vars.scaleValFn{ctrScale,2} = matlabFunction(sym(in.scale{ctrScale,2}),'vars',sym('x',[out.oc.num.states 1]));
	end
end

% Get units and convert them to symbols
unitSym = sym(in.scale(:,1));

% Obtain fieldnames of constants
if isfield(in,'const')
	names = fieldnames(in.const);
else
	names = {};
end

% Determine appropriate scaling
for ctrConst = 1 : 1 : length(names)
	units = sym(in.const.(char(names{ctrConst})){2});
	in.vars.scaleFn.const{ctrConst} = matlabFunction(units,'vars',unitSym);
end

% Obtain fieldnames of constraints
if isfield(in,'constraintVal')
	names = fieldnames(in.constraintVal);
else
	names = {};
end

for ctrConstraint = 1 : 1 : length(names)
	units = sym(in.constraintVal.(char(names{ctrConstraint})){2});
	in.vars.scaleFn.constraintVal{ctrConstraint} = matlabFunction(units,'vars',unitSym);
end

% Scale states, x0, and xf
for ctrState = 1 : 1 : out.oc.num.states
	units = sym(in.oc.state{ctrState,2});
	in.vars.scaleFn.state{ctrState} = matlabFunction(units,'vars',unitSym);
end

% Scale costates (units = COST/STATE)
for ctrCostate = 1 : 1 : out.oc.num.states
	units = sym(['(',in.oc.cost.path{2},')/(',in.oc.state{ctrCostate,2},')']);
	in.vars.scaleFn.costate{ctrCostate} = matlabFunction(units,'vars',unitSym);
end  

% if ~in.convertParametersToStates	
	% Scale independent parameters
	units = in.oc.independentVariable{2};
	in.vars.scaleFn.parameters.independentVariable = matlabFunction(units,'vars',unitSym);
% end

% Scale nu parameters (units = COST/CONSTRAINT)
for ctrConstraint = 1 : 1 : out.oc.num.constraints.initial
	units = ['(',in.oc.cost.path{2},')/(',in.oc.constraint.initial{ctrConstraint,2},')'];
	in.vars.scaleFn.lagrangeMultipliers.initial{ctrConstraint} = matlabFunction(units,'vars',unitSym);
end

for ctrConstraint = 1 : 1 : out.oc.num.constraints.terminal
	units = ['(',in.oc.cost.path{2},')/(',in.oc.constraint.terminal{ctrConstraint,2},')'];
	in.vars.scaleFn.lagrangeMultipliers.terminal{ctrConstraint} = matlabFunction(units,'vars',unitSym);
end

% Scale pii parameters (units = COST/CONSTRAINT))
for ctrPathConstraint = 1 : 1 : out.oc.num.constraints.path
	for ctrN = 1 : 1 : length(out.oc.units.interiorPoint{ctrPathConstraint,1}(:,1))
		
		units = ['(',in.oc.cost.path{2},')/(',char(out.oc.units.interiorPoint{ctrPathConstraint,1}(ctrN,1)),')'];
		in.vars.scaleFn.lagrangeMultipliers.interiorPoint{ctrPathConstraint}{ctrN} = matlabFunction(units,'vars',unitSym);
		
	end
end

% Interior point constraint
for ctrConstraint = 1 : 1 : out.oc.num.constraints.interiorPoint
		
		units = ['(',in.oc.cost.path{2},')/(',char(in.oc.constraint.interiorPoint{ctrConstraint,4}),')'];
		in.vars.scaleFn.lagrangeMultipliers.interiorPoint{out.oc.num.constraints.path+ctrConstraint}{1} = ...
			matlabFunction(units,'vars',unitSym);
end

OUT = in;

end

