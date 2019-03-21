function c = intersect(a,b)
%INTERSECT    Intersection of intervals; empty components set to NaN
%
%   c = intersect(a,b)
%
%Result is an ***outer*** inclusion:
%   c includes the true intersection of a and b
%   if c or components of c are NaN, the (true) intersection of a and b is empty
%
%Input a and b must be both real or both complex
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/06/05     S.M. Rump  correction of result (thanks to Andreas Rauh, Ulm, for hint)
% modified 06/11/06     S.M. Rump  correction of result (thanks to Andreas Rauh, Ulm, for hint)
% modified 09/10/07     S.M. Rump  redesign
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if prod(size(a))>1
    if prod(size(b))==1
      if issparse(a)
        b = b*spones(a);
      else
        b = b*ones(size(a));
      end
    end
  else
    if prod(size(b))>1
      if issparse(b)
        a = a*spones(b);
      else
        a = a*ones(size(b));
      end
    end
  end

  if ~isequal(size(a),size(b))
    error('intersect called with non-matching dimensions')
  end

  if ~isa(a,'intval')
    a = intval(a);
  end
  if ~isa(b,'intval')
    b = intval(b);
  end

  if a.complex & b.complex

    c.complex = 1;
    c.inf = [];
    c.sup = [];
    c.mid = b.mid;
    c.rad = b.rad;

    v = b.mid - intval(a.mid);         % connection of midpoints    
    d2 = sqr(intval(real(v))) + sqr(intval(imag(v)));
    d = sqrt(d2);                      % inclusion of distance of midpoints
    arad = intval(a.rad);
    wng = warning;
    warning off
    sumrad = arad+b.rad;
    index = ( d.inf > sumrad.sup );    % empty intersection
    if any(index(:))
      c.mid(index) = NaN;
      c.rad(index) = NaN;
    end
    x = ( 1 + (arad+b.rad).*(arad-b.rad)./d2 )/2;
    index0 = ( x.inf<=0 );
    if any(index0(:))                  % diameter of a in intersection
      c.mid(index0) = a.mid(index0);
      c.rad(index0) = a.rad(index0);
    end
    index1 = ( x.sup>=1 );
    if any(index1(:))                  % diameter of b in intersection
      c.mid(index1) = b.mid(index1);
      c.rad(index1) = b.rad(index1);
    end
    index = ~( index | index0 | index1 );
    if any(index(:))
      %VVVV  cmid = a.mid(index) + x(index).*v(index);
      %VVVV  c.mid(index) = mid(cmid);
      %VVVV  c.rad(index) = sup( rad(cmid) + sqrt(sqr(arad(index))-sqr(x(index).*d(index))) );
      s.type = '()'; s.subs = {index}; 
      cmid = a.mid(index) + subsref(x,s).*subsref(v,s);
      c.mid = subsasgn(c.mid,s,mid(cmid));
      c.rad = subsasgn(c.rad,s,sup( rad(cmid) + sqrt(sqr(subsref(arad,s))-sqr(subsref(x,s).*subsref(d,s))) ));
      %AAAA  Matlab bug fix      
    end
    warning(wng);
    
    index = ( imag(c.rad)~=0 ) | isnan(a.mid) | isnan(a.rad) | ...
            isnan(b.mid) | isnan(b.rad);
    if any(index(:))
      c.mid(index) = NaN;
      c.rad(index) = NaN;
    end
    
  elseif ~a.complex & ~b.complex    % two real intervals
    
    c.complex = 0;
    c.inf = max(a.inf,b.inf);
    c.sup = min(a.sup,b.sup);
    index = ( c.inf>c.sup ) | isnan(a.inf) | isnan(a.sup) | ...
            isnan(b.inf) | isnan(b.sup);
    if any(index(:))
      c.inf(index) = NaN;
      c.sup(index) = NaN;
    end
    c.mid = [];
    c.rad = [];
    
  else
    error('operands of intersect must be both real or both complex')
  end

  c = class(c,'intval');
 
  setround(rndold)
