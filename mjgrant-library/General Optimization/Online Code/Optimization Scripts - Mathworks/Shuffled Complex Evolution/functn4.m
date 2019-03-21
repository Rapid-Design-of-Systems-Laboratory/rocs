function f=functn(nopt,x)
%
%  This is the Rastrigin Function 
%  Bound: X1=[-1,1], X2=[-1,1]
%  Global Optimum: -2, (0,0)
 
x1 = x(1);
x2 = x(2);
f = x1^2 + x2^2 - cos(18.0*x1) - cos(18.0*x2);

return
