function [x,y] = TwoSum(a,b)
%TWOSUM       Error-free transformation of a+b into x+y with x=fl(a+b)
%
%   [x,y] = TwoSum(a,b)
%
%On return, x+y=a+b and x=fl(a+b) provided no overflow occurs.
%Input a,b may be vectors or matrices as well, single or double precision.
%
%Follows D.E. Knuth: The art of computer programming: seminumerical algorithms, 
%  2nd edition, volume 2, Addison Wesley, 1981.
%Requires 6 flops for scalar input.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  x = a + b;
  if any(~isfinite(x))
    error('overflow occurred in TwoSum')
  end
  z = x - a;
  y = ( a - (x-z) ) + (b-z);
