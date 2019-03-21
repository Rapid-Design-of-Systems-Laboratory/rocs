function f=functn(nopt,x)
%
%  This is the Griewank Function (2-D or 10-D)
%  Bound: X(i)=[-600,600], for i=1,2,...,10
%  Global Optimum: 0, at origin

if nopt==2; d = 200; else d = 4000; end;

u1 = 0.0;
u2 = 1.0;
for j = 1:nopt
    u1 = u1 + x(j)^2 / d;
    u2 = u2 * cos(x(j)/sqrt(j));
end
f = u1 - u2 + 1;

return
