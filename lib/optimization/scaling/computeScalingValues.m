function [scaleVal] = computeScalingValues(in,out,sol)
	% Computes the scaling value for each unit
	
	scaleVal = in.scale;
	for ctr1 = 1 : 1 : length(in.scale(:,1))
		if isnumeric(in.scale{ctr1,2}) % constant scaling factor
			scaleVal{ctr1,2} = in.scale{ctr1,2};
		elseif ~isempty(strfind(in.scale{ctr1,2},'.')) % value from structure in input file
			% Assumes two part structure
			I = strfind(in.scale{ctr1,2},'.');
			struct1 = in.scale{ctr1,2}(1:I-1);
			struct2 = in.scale{ctr1,2}(I+1:end);
			scaleVal{ctr1,2} = in.(struct1).(struct2){1};
		else % comes from trajectory states
			% Determine maximum values of relevant states in scaling
			maxY = num2cell(max(sol.y,[],2));
			% Use previously computed symbolic function
			scaleValFn = in.vars.scaleValFn{ctr1,2};
			scaleVal{ctr1,2} = scaleValFn(maxY{1:out.oc.num.states});
		end
	end	

return

