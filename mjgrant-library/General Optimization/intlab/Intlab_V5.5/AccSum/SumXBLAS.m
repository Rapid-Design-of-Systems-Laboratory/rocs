function [s,t] = SumXBLAS(p)
%SUMXBLAS     Summation 'as if' computed in 2-fold precision
%
%   [s,t] = SumXBLAS(p)
%
%On return, s+t approximates sum(p) with accuracy as if computed 
%  in 2-fold precision with higher order part s and lower order part t.
%
%Implements algorithm BLAS_dsum_x from
%  X. Li, J. Demmel, D. Bailey, G. Henry, Y. Hida, J. Iskandar, 
%    W. Kahan, S. Kang, {S.}, A. Kapur, M. Martin, B. Thompson, {B.},
%    T. Tung, {T.}, D. Yoo: Design, Implementation and Testing of 
%    Extended and Mixed Precision BLAS, ACM Trans. Math. Software, 
%    2(28), p. 152-205, 2002.
%Requires 10n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  s = 0; 
  t = 0;
  for i=1:length(p)
    [t_1,t_2] = TwoSum(s,p(i));
    t_2 = t_2 + t;
    [s,t] = FastTwoSum(t_1,t_2);
  end
