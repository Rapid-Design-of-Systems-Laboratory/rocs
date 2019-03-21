function f=functn(nopt,x)
%
%  This is the Six-hump Camelback Function.
%  Bound: X1=[-5,5], X2=[-5,5]
%  True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
%
x1 = x(1);
x2 = x(2);
f = (4 - 2.1*x1^2 + x1^4/3)*x1^2 + x1*x2 + (-4 + 4*x2^2)*x2^2;
    
return
