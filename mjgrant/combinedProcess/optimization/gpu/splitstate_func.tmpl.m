%codegen
function [ y,xf ] = derivFunc1_split( x0,T0,Tf,h, nSegments, p, constants, constraint, arcSequence)
    assert(isa(constants,'struct'));

	coder.varsize('x0',[__NSTATES__+__MAXARCS__ __MAXARCS__],[1 1]);
    % assert(all(size(x0,1)==[__NSTATES__]));
    assert(isa(x0,'double'));
    assert(isa(T0,'double'));
    assert(isa(Tf,'double'));
    assert(isa(h,'double'));
    assert(isa(nSegments,'double'));

	assert(isa(p,'double'));
	coder.varsize('p',[1000 1],[1 0]);
	coder.varsize('arcSequence',[1 __MAXARCS__]),[0 1];
	assert(isa(arcSequence,'double'));

__CONSTANTS__
    
	nArcs = length(arcSequence);
	tspan = T0:h:Tf;
	nStepsPerArc = ((length(tspan)-1)*4)/nSegments;
	% Allocate memory for all arcs
	y = zeros(nArcs*(nStepsPerArc+1),nSegments*size(x0,1));

	xf = zeros(size(x0,1),nArcs);
	% Repeat for each sub arc
	for arcCtr = 1:nArcs
	    % Get corresponding initial guess and propagate
	    [t,x,xf(:,arcCtr)] = RK4_sub(@__EOMNAME__,tspan,x0(:,arcCtr),(arcCtr),p,constants,constraint, arcSequence);

	    row0 = zeros(1,nSegments*size(x0,1));
	    for k=0:size(x0,1)-1
	        for j=0:nSegments-1
	            row0(k*nSegments+j+1) = x(j*nStepsPerArc+1,k+1);
	        end
	    end
	    % Repeat initial state
	    y((arcCtr-1)*(nStepsPerArc+1)+1,:) = row0;

	    % Copy over remaining states
	    for i=0:nStepsPerArc
	        row = zeros(1,nSegments*size(x0,1));
	        for k=0:size(x0,1)-1
	            for j=0:nSegments-1
	                row(k*nSegments+j+1) = x(nStepsPerArc*j+i+1,k+1);
	            end
	        end
	        y((arcCtr-1)*(nStepsPerArc+1)+i+1,:) = row;
	    end
	end

	%     [t,x,xf] = RK4_sub(@__EOMNAME__,T0:h:Tf,x0,constants,constraint, arcSequence, piiIndex);
	%     nStepsPerArc = ((length(t)-1)*4)/nSegments;
	% y = zeros(nStepsPerArc+1,nSegments*length(x0));
	% 
	% row0 = zeros(1,nSegments*length(x0));
	% for k=0:length(x0)-1
	%     for j=0:nSegments-1
	%         row0(k*nSegments+j+1) = x(j*nStepsPerArc+1,k+1);
	%     end
	% end
	% y(1,:) = row0;
	%     for i=0:nStepsPerArc
	%     row = zeros(1,nSegments*length(x0));
	%     for k=0:length(x0)-1
	%         for j=0:nSegments-1
	%             row(k*nSegments+j+1) = x(nStepsPerArc*j+i+1,k+1);
	%         end
	%     end
	%     y(i+1,:) = row;
	%     end
end

