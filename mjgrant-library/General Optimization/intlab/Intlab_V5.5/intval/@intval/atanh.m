function y = atanh(x)
%ATANH        Implements  atanh(x)  for intervals
%
%   y = atanh(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, sparse input,
%                                  major revision, improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  'realstdfctsexcptnignore' added and some
%                                     improvements, tocmplx replaced by cintval
% modified 09/06/07     S.M. Rump  approximate std fcts removed, exceptional arguments
%

  if issparse(x)
    [ix,jx,sx] = find(x);
    [m,n] = size(x);
    y = sparse(ix,jx,atanh(full(sx)),m,n);
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
    y = log( 2./(1-x) - 1 )/2;
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
    index = ( x.sup<-1 ) | ( x.inf>1 );   % completely exceptional indices
  else
    index = ( x.inf < -1 ) | ( x.sup > 1 ) ;
    if ( INTLAB_INTVAL_STDFCTS_EXCPTN<2 ) & any(index(:))
      if INTLAB_INTVAL_STDFCTS_EXCPTN==1
        warning('ATANH: Real interval input out of range changed to be complex')
      end
      y = x;
      %VVVV  y(index) = atanh(cintval(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,atanh(cintval(subsref(x,s))));
      %AAAA  Matlab bug fix
      index = ~index;
      if any(index(:))
        %VVVV  y(index) = atanh(x(index));
        s.type = '()'; s.subs = {index}; y = subsasgn(y,s,atanh(subsref(x,s)));
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
  wng = warning;
  warning off

  % treat non-exceptional arguments
  xinf = x.inf(:);
  xsup = x.sup(:);

  IndexInfPos = ( xinf>=0 );
  len1 = sum(IndexInfPos);
  IndexSupNeg = ( xsup<=0 );
  len2 = sum(IndexSupNeg);

  Y = atanh_pos( [ xinf(IndexInfPos) ; -xsup(IndexSupNeg) ] , -1 );
  y.inf(IndexInfPos) = Y(1:len1);
  y.sup(IndexSupNeg) = -Y( len1+1 : end );

  IndexInfNeg = ( xinf<0 );
  len1 = sum(IndexInfNeg);
  IndexSupPos = ( xsup>0 );
  len2 = sum(IndexSupPos);

  Y = atanh_pos( [ -xinf(IndexInfNeg) ; xsup(IndexSupPos) ] , 1 );
  y.inf(IndexInfNeg) = -Y(1:len1);
  y.sup(IndexSupPos) = Y( len1+1 : len1+len2 );

  y.inf(xinf==-1) = -inf;
  y.sup(xsup==-1) = -inf;
  y.inf(xinf==1) = inf;
  y.sup(xsup==1) = inf;

  if any(index(:))                       % exceptional arguments
    y.inf(index) = NaN;
    y.sup(index) = NaN;
  end

  setround(rndold)
  warning(wng)
  