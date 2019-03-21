function res = Dot2(x,y)
%DOT2         Dot product 'as if' computed in 2-fold (quadruple) precision
%
%   res = Dot2(x,y)
%
%On return, res approximates x'*y with accuracy as if computed 
%  in 2-fold precision.
%
%Implements algorithm Dot2 from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires 25n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  [p,s] = TwoProduct(x(1),y(1));
  for i=2:length(x)
    [h,r] = TwoProduct(x(i),y(i));
    [p,q] = TwoSum(p,h);
    s = s + ( q + r );
  end
  res = p + s;
  