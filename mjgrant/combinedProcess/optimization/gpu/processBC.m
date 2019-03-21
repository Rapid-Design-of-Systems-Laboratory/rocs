function [M, N] = processBC(bcfunc,in,constants,varargin)
	
	YL = sym('X0',[2*in.oc.num.states+1 1]);
	YR = sym('XF',[2*in.oc.num.states+1 1]);
	
	numNu0 = in.oc.num.constraints.initial;
	numNuF = in.oc.num.constraints.terminal;
	
	p = sym('p',[numNu0+numNuF 1]);
    % Create symbolic variables for the auxiliary data

	if isfield(in,'constraintVal')
		constraints = in.constraintVal;
	else
		constraints = struct();
	end
	res = bcfunc(YL,YR,p,constants,constraints,varargin{:});
	
	M = double(jacobian(res,[YL;p]));
	N = double(jacobian(res,[YR;sym('dummy_',[numNu0+numNuF,1])]));
	% P = double(jacobian(res,p))
	% keyboard;
end