function res = DotK(x,y,K)
%DOTK         Dot product 'as if' computed in K-fold (quadruple) precision
%
%   res = DotK(x,y,K)
%
%On return, res approximates x'*y with accuracy as if computed 
%  in K-fold precision.
%
%Implements algorithm DotK from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires (12K+1)n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  n = length(x);
  r = zeros(2*n,1);
  [p,r(1)] = TwoProduct(x(1),y(1));
  for i=2:n
    [h,r(i)] = TwoProduct(x(i),y(i));
    [p,r(n+i-1)] = TwoSum(p,h);
  end
  r(2*n) = p;
  res = SumK(r,K-1);
  