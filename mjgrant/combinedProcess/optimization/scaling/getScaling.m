function [scale] = getScaling(units,scaleVal)
	% This function determine the appropriate scaling for a particular parameter.
	
	scale = sym(units);

	% substitute units one by one
	for ctr2 = 1 : 1 : length(scaleVal(:,1))
		scale = subs(scale,sym(scaleVal{ctr2,1}),sym(scaleVal{ctr2,2}))
	end

	scale = double(scale);
return