function u = rdivide(a,b)
%RDIVIDE      Slope elementwise right division  a ./ b
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  global INTLAB_SLOPE

  if ~isa(a,'slope')
    a = slope(a);
  end
  if ~isa(b,'slope')
    b = slope(b);
  end

  na = size(a.r.inf,1);
  nb = size(b.r.inf,1);

  if ( na==1 ) & ( nb~=1 )
    a.r = a.r(ones(nb,1),:);
    a.s = a.s(ones(nb,1),:);
    u.size = b.size;
  elseif ( na~=1 ) & ( nb==1 )
    b.r = b.r(ones(na,1),:);
    b.s = b.s(ones(na,1),:);
    u.size = a.size;
  else
    if ~isequal(a.size,b.size)
      error('dimensions not compatible for minus')
    end
    u.size = a.size;
  end

  u.r = a.r ./ b.r;
  indexc = 1:INTLAB_SLOPE.NUMVAR;
  indexr = 2:INTLAB_SLOPE.NUMVAR+1;
  u.s = intersect( ( a.s - u.r(:,indexc) .* b.s ) ./ b.r(:,indexr) , ...
                   ( a.s - u.r(:,indexr) .* b.s ) ./ b.r(:,indexc) );

  u.r = rangeimprove(u);

  u = class(u,'slope');
  
  if rndold~=0
    setround(rndold)
  end
