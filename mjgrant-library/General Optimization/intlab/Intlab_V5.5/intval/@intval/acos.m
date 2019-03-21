function y = acos(x)
%ACOS         Implements  acos(x)  for intervals
%
%   y = acos(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, sparse input,
%                                  pos/neg split, major revision,
%                                  improved accuracy, corrected
%                                  branchcut
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  'realstdfctsexcptnignore' added and
%                                     some improvements
% modified 09/06/07     S.M. Rump  approximate std fcts removed, exceptional arguments
%

  if issparse(x)
    index = ( x==0 );
    if any(index(:))                    % treat zero indices
      global INTLAB_INTVAL_STDFCTS_PI
      PI2 = intval(INTLAB_INTVAL_STDFCTS_PI.PI2INF,INTLAB_INTVAL_STDFCTS_PI.PI2SUP,'infsup');
      y = intval(repmat(PI2,size(x)));
      index = ~index;
      %VVVV  y(index) = acos(full(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,acos(full(subsref(x,s))));
      %AAAA  Matlab bug fix
    else
      y = acos(full(x));
    end
    return
  end
      
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if x.complex
%   y = -i * log( x + sqrt(x.^2-1) );
    y = i * acosh(x);
    index = ( real(y.mid) < 0 );
    y.mid(index) = -y.mid(index);
    if rndold~=0
      setround(rndold)
    end
    return
  end

  % input x real and full
  % real range of definition:  [-1,1]
  global INTLAB_INTVAL_STDFCTS_EXCPTN   
  if INTLAB_INTVAL_STDFCTS_EXCPTN==3    % ignore input out of range
    global INTLAB_INTVAL_STDFCTS_EXCPTN_   
    index = ( x.inf<-1 );               % (partially) exceptional indices
    if any(index(:))
      x.inf(index) = -1;
      INTLAB_INTVAL_STDFCTS_EXCPTN_ = 1;
    end
    index = ( x.sup>1 );                % (partially) exceptional indices
    if any(index(:))
      x.sup(index) = 1;
      INTLAB_INTVAL_STDFCTS_EXCPTN_ = 1;
    end
    exceptions = ( x.sup < -1 ) | ( x.inf > 1 );   % completely exceptional indices
  else
    exceptions = ( x.inf < -1 ) | ( x.sup > 1 );
    if ( INTLAB_INTVAL_STDFCTS_EXCPTN<2 ) & any(exceptions(:))
      if INTLAB_INTVAL_STDFCTS_EXCPTN==1
        warning('ACOS: Real interval input out of range changed to be complex')
      end
      y = x;
      %VVVV  y(exceptions) = acos(cintval(x(exceptions)));
      s.type = '()'; s.subs = {exceptions}; y = subsasgn(y,s,acos(cintval(subsref(x,s))));
      %AAAA  Matlab bug fix
      index = ~exceptions;
      if any(index(:))
        %VVVV  y(index) = acos(x(index));
        s.type = '()'; s.subs = {index}; y = subsasgn(y,s,acos(subsref(x,s)));
        %AAAA  Matlab bug fix
      end
      if rndold~=0
        setround(rndold)
      end
      return
    end
  end

  % input x real and full
  y = x;

  % treat non-exceptional cases
  xinf = x.inf(:);
  xsup = x.sup(:);

  % switch off warning since exceptional values may occur
  wng = warning;
  warning off

  index = ( xsup>=0 );
  if any(index)
    y.inf(index) = - asin_pos_( xsup(index) , 1 , -1 );
  end

  index = ( xsup<0 );
  if any(index)
    y.inf(index) = asin_pos_( -xsup(index) , -1 , 1 );
  end

  index = ( xinf>=0 );
  if any(index)
    y.sup(index) = - asin_pos_( xinf(index) , -1 , -1 );
  end

  index = ( xinf<0 );
  if any(index)
    y.sup(index) = asin_pos_( -xinf(index) , 1 , 1 );
  end

  if any(exceptions(:))                   % exceptional arguments
    y.inf(exceptions) = NaN;
    y.sup(exceptions) = NaN;
  end

  % restore warning status
  warning(wng);

  y.inf = reshape(y.inf,size(x.inf));
  y.sup = reshape(y.sup,size(x.sup));

  setround(rndold)
