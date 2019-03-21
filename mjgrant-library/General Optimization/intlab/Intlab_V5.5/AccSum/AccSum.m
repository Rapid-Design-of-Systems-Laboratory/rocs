function res = AccSum(p)
%ACCSUM       Faithful rounding of sum(p)
%
%   res = AccSum(p)
%
%On return, res is a faithful rounding of sum(p), also in the presence
%  of underflow. Input vector p may be single or double precision.
%
%Implements Algorithm 4.5 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
%    Faithful Rounding, to appear. 
%Requires (4m+3)n flops for m executions of repeat-until loop in Transform.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  if isa(p,'double')
    nmax = 2^26-2;          % nmax = 67,108,864
  else
    nmax = 2^12-2;          % nmax = 4,094
  end
  if length(p)>nmax
    error(['maximum length of input vector ' int2str(nmax) '; use AccSumHugeN'])
  end
  
  [tau1,tau2,p] = Transform(p,0);
  res = tau1 + ( tau2 + sum(p) );
