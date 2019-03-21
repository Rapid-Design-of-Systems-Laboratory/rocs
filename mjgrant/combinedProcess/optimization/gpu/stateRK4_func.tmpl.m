%codegen
function [t,x] = derivFunc1_RK4(x0,nTimeStepsPerArc,p, constants, constraint, arcSequence)
	assert(isa(constants,'struct'));
	assert(isa(x0,'double'));
	coder.varsize('x0',[__NSTATES__+__MAXARCS__ __MAXARCS__],[1 1]);
	assert(isa(nTimeStepsPerArc,'double'));
	
	assert(isa(p,'double'));
	coder.varsize('p',[1000 1],[1 0]);
	coder.varsize('arcSequence',[1 __MAXARCS__]),[0 1];
	assert(isa(arcSequence,'double'));

__CONSTANTS__
	h = 1.0/nTimeStepsPerArc;

	x = NaN(size(x0,1),(nTimeStepsPerArc+1)*length(arcSequence));
	t = NaN((nTimeStepsPerArc+1)*length(arcSequence),1);

	T0 = 0;
	Tf = 1;
	for arcCtr = 1 : 1 : length(arcSequence)
		tspan = T0:h:Tf;
		[t((arcCtr-1)*length(tspan)+1:arcCtr*length(tspan),1),x(:,(arcCtr-1)*length(tspan)+1:arcCtr*length(tspan))] = ...
					RK4(@__EOMNAME__,tspan,x0(:,arcCtr),arcCtr,p,constants,constraint, arcSequence);
		T0 = T0 + 1;
		Tf = Tf + 1;
	end
end

