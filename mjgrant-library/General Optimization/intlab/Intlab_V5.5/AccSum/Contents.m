%Reference implementations for accurate summation and dot product algorithms
%
%New summation, dot product and related routines: high accuracy
%  AccSum        - faithful rounding of sum(p_i)
%  AccSumHugeN   - faithful rounding of sum(p_i) for large dimension
%  AccSumK       - K-fold faithful rounding of sum(p_i)
%  AccSign       - sign of sum(p_i)
%  NearSum       - sum(p_i) rounded to nearest
%  DownSum       - sum(p_i) rounded downwards
%  UpSum         - sum(p_i) rounded upwards
%
%
%New summation, dot product and related routines: high precision
%  Sum2          - summation with quad precision
%  SumK          - summation with K-fold precision
%  Dot2          - dot product with quad precision
%  Dot2Err       - dot2 with rigorous error bounds without directed rounding
%  DotK          - dot product with K-fold precision
%
%
%New and fast summation algorithms with high accuracy
%  PrecSum       - accurate and fast up to large condition number
%
%
%Error-free transformations
%  TwoSum        - transformation of a+b into x+y with x=fl(a+b)
%  FastTwoSum    - as TwoSum provided input is ordered in absolute value
%  TwoProduct    - transformation of a*b into x+y with x=fl(a*b)
%  Split         - transformation of a into two 'halves' x+y
%  ExtractVector - extract higher and lower part of vector
%  Transform     - Transformation of vector into high part and low order vector
%
%
%Utility routines for new summation and dot product routines
%  NextPowerTwo  - next power of 2 of integer without branch
%
%
%Reference implementations of competitors
%  SumXBLAS      - summation as in XBLAS
%  DotXBLAS      - dot product as in XBLAS
%  PriestSum     - Priest's doubly compensated summation
%  ZDSum         - Zielke/Drygalla summation
%
%

% written  05/07/08     S.M. Rump
%

%
%New algorithms based on
%
%Accurate summation and dot product with specified precision
%  T. Ogita, S.M. Rump, and S. Oishi. Accurate Sum and Dot Product, 
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005.
%
%Accurate summation and dot product with specified accuracy
%  S.M. Rump, T. Ogita, and S. Oishi. Accurate Floating-point Summation I: 
%     Faithful Rounding. accepted for publication in SISC, 2005-2008. 
%  S.M. Rump, T. Ogita, and S. Oishi. Accurate Floating-point Summation II: 
%     Sign, K-fold Faithful and Rounding to Nearest. submitted for publication 
%     in SISC, 2005-2008. 
%
%Fast and high precision summation and dot products
%  S.M. Rump, T. Ogita, and S. Oishi. Fast High Precision Summation, 
%     submitted for publication.

%
% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
