function res = SumK(p,K)
%SUMK         Summation 'as if' computed in K-fold precision
%
%   res = SumK(p,K)
%
%On return, res approximates sum(p) with accuracy as if computed 
%  in K-fold precision.
%
%Implements algorithm SumK from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires (6K-5)n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  for k=1:K-1
    for i=2:length(p)
      [p(i),p(i-1)] = TwoSum(p(i),p(i-1));
    end
  end
  res = sum(p(1:end-1)) + p(end);
