function [out_oc] = addExtraStates(oc,numArcs,issym)
	out_oc = oc;
	if nargin < 3
		issym = false;
	end
	if ~issym
		numStates = length(oc.state(:,1));
	else
		numStates = length(oc.state.var);
	end

	for ctr = 1 : 1 : numArcs
		if ~issym
			out_oc.state{numStates+ctr,1} = strcat(['t',int2str(ctr)]);
			out_oc.state{numStates+ctr,2}= 's';
			out_oc.stateRate{numStates+ctr,1} = '0';
		else
			out_oc.state.var(numStates+ctr) = sym(['t',int2str(ctr)]);
			out_oc.state.rate(numStates+ctr) = 0;
		end
	end
end