function [sol] = newtonRaphson(func,x0,tol,varargin)
    assert(isa(x0,'double'));
    assert(isa(func,'function_handle'));
    assert(isa(tol,'double'));
    % coder.varsize('x0',[100 1],[1 0]);
    
    N = length(x0);
    
    x = x0;
    fx = inf(N,1);  % Function value
    while true
        fx = func(x,varargin{:});
        if max(abs(fx)) <= tol
            break;
        end
        J  = finiteDifferenceJacobian(func,x0,varargin{:});
        x1 = x - J\fx;
             
        x = x1;
    end
    sol = x;
end