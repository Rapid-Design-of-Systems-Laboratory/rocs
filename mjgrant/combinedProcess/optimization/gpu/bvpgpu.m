function [sol,jac] = bvpgpu(eomfunc,bcfunc,solInit,in,out,constants,constraint,arcSequence,varargin)

	% Assume files already created by runCombinedProcess

	% Copy initial guess
	y = solInit.y;

	repVals = find(diff(solInit.x) == 0);
	repVals = [0 repVals];
	
	if length(arcSequence) > 3
		% in.gpuSolve.TimeSteps = 512*3;
		% in.gpuSolve.dampingFactor = 0.8;
	end
	for i=1:length(arcSequence)
		x0g(:,i) = y(:,repVals(i)+1);
	end
	oldx = solInit.x;
	% solInit.x = tspan;
% 
% 	repVals = [repVals length(solInit.x)];
% 	
% 	if length(repVals) ~= length(arcSequence)
% 		error('Number of arcs in input does not match actual number of arcs');
% 	end
% 	for i=1:length(arcSequence)
% 		for j=1:in.gpuSolve.numSubArcs
% 			if i == 1
% 				subArcPos = ceil(j*repVals(i)/in.gpuSolve.numSubArcs);
% 			else
% 				subArcPos = repVals(i-1)+ ceil(j*(repVals(i) - repVals(i-1))/in.gpuSolve.numSubArcs);
% 			end
% 			x0g(:,(i-1)*in.gpuSolve.numSubArcs + j) = solInit.y(:,subArcPos);
% 
% 		end
% 	end
	
	t0 = solInit.x(1);
	% tf = solInit.x(end);

	constArray = cell2mat(struct2cell(constants));
	i = 0;

	paramGuess = solInit.parameters;
	numArcs = length(arcSequence);
	% zeroPadding = zeros((in.oc.num.states*2+numArcs),in.oc.num.parameters);
	% stm00 = eye(14,14);
  %warning('error', 'MATLAB:nearlySingularMatrix');
  % warning('error', 'MATLAB:singularMatrix');
	% last_dx0 = inf(size(x0g,1)*numArcs+length(paramGuess),1);
	while true
		[phi,xf] = propagate_stm(eomfunc,t0,1.0,in.gpuSolve.TimeSteps,x0g,paramGuess,constants,constArray,in,constraint,arcSequence,varargin{:});
		format short;
    % xf
    if length(arcSequence) > 3
      % keyboard;
    end
        res = bcfunc(x0g,xf,paramGuess,constants,constraint,arcSequence,varargin{:});
    % res
% keyboard;
        % Augmented Jacobian
		jac = bcJac(bcfunc,phi,in,x0g,xf,paramGuess,constants,constraint,arcSequence,varargin{:});
		% if length(arcSequence) > 2
% %             keyboard
            % max(abs(res))
        % end
		if any(isnan(res))
			error('NaN in residue vector. Shooting method diverged. Stopping.');
		end
		if max(abs(res)) < in.gpuSolve.Tol
			break;
		end
		% MATLAB:nearlySingularMatrix, MATLAB:singularMatrix
		try
			dx0 = -jac\res;
		catch
			dx0 = pinv(jac)*-res;
		end
		% while max(abs(dx0)) > max(abs(last_dx0))
		% 	dx0 = dx0*0.9;
		% end
		dxIndex = 1;
		for arcCtr = 1:length(arcSequence)
			x0g(:,arcCtr) = x0g(:,arcCtr) + in.gpuSolve.dampingFactor*dx0(dxIndex:dxIndex+size(x0g,1)-1);
			dxIndex = dxIndex + size(x0g,1);
		end
		paramGuess = paramGuess	+ in.gpuSolve.dampingFactor*dx0(dxIndex:dxIndex+length(paramGuess)-1);
		% [res, dx0]
		% max(abs(res))
		if i>in.gpuSolve.maxIterations
			error('Maximum iterations exceeded!');
			break;
		end
		i = i+1;
	end	
	[sol.x,sol.y] = derivFuncRegion_RK4_mex(x0g,in.gpuSolve.TimeSteps,[1],constants,constraint,arcSequence);
	sol.x = sol.x';
	sol.parameters = paramGuess;
end
