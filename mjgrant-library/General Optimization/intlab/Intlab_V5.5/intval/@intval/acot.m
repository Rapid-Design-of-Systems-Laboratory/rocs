function y = acot(x)
%ACOT         Implements  acot(x)  for intervals
%
%   y = acot(x)
%
%interval standard function implementation
%

% written  10/16/98     S.M. Rump
% modified 06/24/99     S.M. Rump  complex allowed, sparse input,
%                                  major revision, improved accuracy
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/06/07     S.M. Rump  improved performance
%

  if issparse(x)
    index = ( x==0 );
    if any(index(:))                    % treat zero indices
      global INTLAB_INTVAL_STDFCTS_PI
      PI2 = intval(INTLAB_INTVAL_STDFCTS_PI.PI2INF,INTLAB_INTVAL_STDFCTS_PI.PI2SUP,'infsup');
      y = intval(repmat(PI2,size(x)));
      index = ~index;
      %VVVV  y(index) = acot(full(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,acot(full(subsref(x,s))));
      %AAAA  Matlab bug fix
    else
      y = acot(full(x));
    end
    return
  end
      
  y = atan(1./x);
