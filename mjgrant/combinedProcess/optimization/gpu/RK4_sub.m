%codegen
function [tout,xout,xf] = RK4_sub(f,tvect,x0,varargin)
    assert(isa(f,'function_handle'));
    assert(isa(x0,'double'));
    assert(isa(tvect,'double'));
    
    n=1;
    t = tvect;
    t(1)=tvect(1);
    [a b]=size(tvect);

    if a>b
        n_max=a;
    else
        n_max=b;
    end
    h = zeros(size(tvect));
    x = zeros(length(x0),n_max);
    x(:,1)=x0;
    xout = zeros(length(x0),(n_max-1)*4+1);
    xout(:,1) = x0;
    
    k = zeros(length(x0),4);
    for n=1:n_max-1

            h(n)=tvect(n+1)-tvect(n);

            k(:,1)=h(n)*f(t(n),x(:,n),varargin{:});
            k(:,2)=h(n)*f(t(n)+h(n)/2,x(:,n)+k(:,1)/2,varargin{:});
            k(:,3)=h(n)*f(t(n)+h(n)/2,x(:,n)+k(:,2)/2,varargin{:});
            k(:,4)=h(n)*f(t(n)+h(n),x(:,n)+k(:,3),varargin{:});
            x(:,n+1)=x(:,n)+(k(:,1)+2*k(:,2)+2*k(:,3)+k(:,4))/6;

            xout(:,(n-1)*4+2) = x(:,n)+k(:,1)/2;
            xout(:,(n-1)*4+3) = x(:,n)+k(:,2)/2;
            xout(:,(n-1)*4+4) = x(:,n)+k(:,3);
            xout(:,(n*4+1)) = x(:,n+1);
            t(n+1)=t(n)+h(n);
            n1=n;

    end
    tout=t';
    xout = xout';
	xf = x(:,end);
    %xout=x';
end
% function [tout,xout] = RK4_sub(f,tvect,x0,varargin)
% 
% n=1;
% t(1)=tvect(1);
% [a b]=size(tvect);
% 
% if a>b
%     n_max=a;
% else
%     n_max=b;
% end
% 
% x(:,1)=x0;
% xout = zeros(length(x0),(n_max-1)*4+1);
% xout(:,1) = x0;
% for n=1:n_max-1
%       
%         h(n)=tvect(n+1)-tvect(n);
%         
%         k(:,1)=h(n)*f(t(n),x(:,n),varargin{:});
%         k(:,2)=h(n)*f(t(n)+h(n)/2,x(:,n)+k(:,1)/2,varargin{:});
%         k(:,3)=h(n)*f(t(n)+h(n)/2,x(:,n)+k(:,2)/2,varargin{:});
%         k(:,4)=h(n)*f(t(n)+h(n),x(:,n)+k(:,3),varargin{:});
%         x(:,n+1)=x(:,n)+(k(:,1)+2*k(:,2)+2*k(:,3)+k(:,4))/6;
%         
% 		xout(:,(n-1)*4+2) = x(:,n)+k(:,1)/2;
% 		xout(:,(n-1)*4+3) = x(:,n)+k(:,2)/2;
% 		xout(:,(n-1)*4+4) = x(:,n)+k(:,3);
% 		xout(:,(n*4+1)) = x(:,n+1);
%         t(n+1)=t(n)+h(n);
%         n1=n;
%     
% end
% tout=t';
% xout = xout';
% %xout=x';
% 
% end