function [J] = finiteDifferenceJacobian(func,x0,varargin)
    assert(isa(func,'function_handle'));
    assert(isa(x0,'double'));
    % coder.varsize('x0',[100 1],[1 0]);
    
    N = length(x0);
    h = 1e-6;
    J = zeros(N,N);
    x = x0;
    for i=1:N
		x(i,1) = x(i,1) + h;
		fxh = func(x,varargin{:});
		x(i,1) = x(i,1) - h;
		dx = fxh/h - func(x,varargin{:})/h;
		J(:,i) = (dx);        
    end
end