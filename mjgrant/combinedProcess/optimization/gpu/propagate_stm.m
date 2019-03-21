function [phi,xf] = propagate_stm(eomfunc,T0,Tf,nStepsPerArc,x0,p,constants,constArray,in,constraint,arcSequence,varargin)

	config = in.gpuSolve;
	
	if config.sttOrder ~= 1
		err = MException('PropagateSTM:sttOrder', 'propagate_stm() works only when config.sttOrder == 1');
		throw(err);
		return;
	end
	
	h  = (Tf-T0)/nStepsPerArc;
	
	nSegments = min(config.ArcSegments,nStepsPerArc);
	allowed_segments = [1,2,4,8,16,32,48,64,128,256,512];
	if isfield(config,'ArcSegments') && sum(find(allowed_segments == config.ArcSegments))>0
		nSegments = min(config.ArcSegments,config.TimeSteps);
	else
		nSegments = 32;
	end
	
	
	constraint_fields = fieldnames(constraint);
	constraintArray = NaN(length(constraint_fields),1);
	for ctr=1:length(constraint_fields)
		fieldname = char(constraint_fields(ctr));
		constraintArray(ctr) = constraint.(fieldname);
	end

	% newArcSequence = NaN(1,length(arcSequence)*2);
	% for i=1:length(arcSequence)
	% 	for j=1:in.gpuSolve.numSubArcs
	% 		newArcSequence(in.gpuSolve.numSubArcs*i-(j-1)) = arcSequence(i);
	% 	end
	% end

	% keyboard
	% EOM files already exist
	% eomname = func2str(eomfunc);
% constraintArray
	[y,xf] = derivFuncRegion_split_mex(x0,T0,Tf,h,nSegments,p,constants,constraint,arcSequence);

	% [y,xf] = derivFuncRegion_split_mex(x0,T0,Tf,h,nSegments,p,constants,constraint,newArcSequence);
	% toc
	eomfile = config.eomfile;
	% phi = propagator_stm_gpu(eomfile,y',T0,Tf,nSteps,[1],constArray,in.oc.num.states,nSegments);
	% tic

	phi = propagator_multi(eomfile,y',T0,Tf,nStepsPerArc,constArray,constraintArray,in.oc.num.states*2+length(arcSequence),arcSequence,nSegments);
	% toc
	% stmprop_time = toc;
	
	% fprintf('STM propagation on GPU : %0.4f s\n',stmprop_time);
end