function c = mid(a)
%MID          Implements  mid(a)  for intervals (rounded)
%
%   c = mid(a)
%
% mid(a) and rad(a) computed such that
%    alpha  in  < mid(a) , rad(a) >  for all alpha in a
%

% written  10/16/98     S.M. Rump
% modified 06/22/99     S.M. Rump  for sparse matrices
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 09/10/07     S.M. Rump  performance, huge arrays
%

  if a.complex
    c = a.mid;
  else
    e = 1e-30;
    if 1+e==1-e                           % fast check for rounding to nearest
      rndold = 0;
    else
      rndold = getround;
    end
    setround(1)
    c = 0.5*a.inf+0.5*a.sup;
    setround(rndold)                    % set rounding to previous value
    index = isinf(a.inf) & isinf(a.sup) & isnan(c);
    [m,n] = size(a.inf);
    if m*n<2^31
      if any(index(:))
        c(index) = 0;
      end
    else            % take care of huge matrices
      index_ = any(index);
      if any(index_(:))
        [I1,J1,ainf] = find(a.inf);
        [I2,J2,asup] = find(a.sup);
        [I,J,ainfsup] = find(sparse([I1(:);I2(:)],[J1(:);J2(:)],[ainf(:);complex(0,asup(:))],m,n));
        index = isinf(real(ainfsup)) & isinf(imag(ainfsup)) & isnan(real(ainfsup)+imag(ainfsup));
        ainfsup(index) = 0;
        setround(1)
        c = sparse(I,J,0.5*real(ainfsup)+0.5*imag(ainfsup),m,n);
        setround(rndold)                % set rounding to previous value
      end
    end
  end
