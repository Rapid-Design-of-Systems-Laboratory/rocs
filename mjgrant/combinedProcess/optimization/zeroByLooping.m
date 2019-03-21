function variableSave2 = zeroByLooping(functionHandle,interval,varargin)
	coder.extrinsic('keyboard');

	variableSave = NaN;
	zeroValSave = inf;
	
	% interval(1)
	% interval(2)
	for variable = linspace(interval(1),interval(2),100)
		
		zeroVal = alfaSolutionFunction(variable,varargin{:});
		if abs(zeroVal) < zeroValSave
			zeroValSave = abs(zeroVal);
			variableSave = variable;
		end
		
	end
	
	if isnan(variableSave) || isinf(variableSave)
		
		varargin{1}
		% varargin{4}
		variableSave
		
	end
	
	[variableSave2,zeroValSave2] = fzero(@alfaSolutionFunction,variableSave,[],varargin{:});
	
	% zeroValSave = zeroValSave
	% variableSave = variableSave
	% 
	% zeroValSave2 = zeroValSave2
	% variableSave2 = variableSave2
	
	% keyboard

return

