function f=functn(nopt,x)
%
%  This is the Rosenbrock Function
%  Bound: X1=[-5,5], X2=[-2,8]
%  Global Optimum: 0,(1,1)
%
x1 = x(1);
x2 = x(2);
a = 100;
f = a * (x2 - x1^2)^2 + (1 - x1)^2;

return
 