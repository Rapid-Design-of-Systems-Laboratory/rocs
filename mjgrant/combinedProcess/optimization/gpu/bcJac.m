function [jac] = bcJac(bcFunc,phiList,in,YL,YR,p,const,constraint,arcSequence,varargin)
	% phiList has phi's in order of arcs, [phi1,phi2,...]

	nArcs = length(arcSequence);	
	nStatesAndCostates = 2*in.oc.num.states;
	nTotalODEs = nStatesAndCostates + nArcs;
	nParams = length(p);

	jac = NaN(nTotalODEs*nArcs+length(p));
	if strcmpi(in.gpuSolve.derivativeMethod,'csd')
		h = 1e-50;
	
		YL_c = complex(YL);
		YR_c = complex(YR);
		p_c  = complex(p);

		jacCol = 1;
		for arcIndex = 1:nArcs
			phi = phiList(:,(arcIndex-1)*nTotalODEs+1:arcIndex*nTotalODEs);
			M = NaN(nTotalODEs*nArcs+nParams,nTotalODEs);
			N = NaN(nTotalODEs*nArcs+nParams,nTotalODEs);
			for ctr=1:nTotalODEs
				YL_c(ctr,arcIndex) = YL_c(ctr,arcIndex) + 1i*h;
				res_c = bc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:});
				YL_c(ctr,arcIndex) = YL_c(ctr,arcIndex) - 1i*h;	
				M(:,ctr) = imag(res_c)/h;
			end
			for ctr=1:nTotalODEs
				YR_c(ctr,arcIndex) = YR_c(ctr,arcIndex) + 1i*h;
				res_c = bc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:});
				YR_c(ctr,arcIndex) = YR_c(ctr,arcIndex) - 1i*h;
				N(:,ctr) = imag(res_c)/h;
			end
			arcJacobian = M+N*phi;
			jac(:,jacCol:jacCol+nTotalODEs-1) = arcJacobian;
			jacCol = jacCol + nTotalODEs;
		end
		for ctr=1:nParams
			p_c(ctr) = p_c(ctr) + 1i*h;
			res_c = bc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:});
			p_c(ctr) = p_c(ctr) - 1i*h;
			jac(:,jacCol) = imag(res_c)/h;
			jacCol = jacCol + 1;
		end
	else
		h = 1e-8;
		YL_c = (YL);
		YR_c = (YR);
		p_c  = (p);
	
		jacCol = 1;
		for arcIndex = 1:nArcs
			phi = phiList(:,(arcIndex-1)*nTotalODEs+1:arcIndex*nTotalODEs);
			M = NaN(nTotalODEs*nArcs+nParams,nTotalODEs);
			N = NaN(nTotalODEs*nArcs+nParams,nTotalODEs);
			for ctr=1:nTotalODEs
				YL_c(ctr,arcIndex) = YL_c(ctr,arcIndex) + h;
				res = bcFunc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:});
				YL_c(ctr,arcIndex) = YL_c(ctr,arcIndex) - h;		
				res = res/h - bcFunc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:})/h;
				M(:,ctr) = (res);
			end
			for ctr=1:nTotalODEs
				YR_c(ctr,arcIndex) = YR_c(ctr,arcIndex) + h;
				res = bcFunc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:});
				YR_c(ctr,arcIndex) = YR_c(ctr,arcIndex) - h;		
				res = res/h - bcFunc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:})/h;
				N(:,ctr) = (res);
			end
			arcJacobian = M+N*phi;
			jac(:,jacCol:jacCol+nTotalODEs-1) = arcJacobian;
			jacCol = jacCol + nTotalODEs;
		end
		for ctr=1:nParams
			p_c(ctr) = p_c(ctr) + h;
			res = bcFunc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:});
			p_c(ctr) = p_c(ctr) - h;
			res = res/h - bcFunc(YL_c,YR_c,p_c,const,constraint,arcSequence,varargin{:})/h;
			jac(:,jacCol) = res;
			jacCol = jacCol + 1;
		end
		% rank(jac)
	end
return;
end