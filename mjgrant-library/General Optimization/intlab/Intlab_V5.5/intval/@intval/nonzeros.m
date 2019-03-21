function x = nonzeros(a)
%NONZEROS     Implements  nonzeros(a)  for sparse interval matrix
%
%   x = nonzeros(a)
%
%Functionality as in Matlab.
%

% written  08/09/02     S.M. Rump 
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/25/07     S.M. Rump  huge indices for sparse matrices
% modified 02/18/08     S.M. Rump  no lolumn vector output for ND arrays
%

  if a.complex
    if isequal(a.rad,0)
      x = cintval(nonzeros(a.mid));
    else
      I = find(spones(a.mid)+spones(a.rad));
      x = intval(full(a.mid(I)),full(a.rad(I)),'midrad');
      if size(x.mid,1)==1
        x.mid = x.mid.';
        x.rad = x.rad';
      end
    end
  else
    if issparse(a.inf)              % take care of huge indices
      [I,J,asup] = find(a.sup);
      [m,n] = size(a.sup);
      x = nonzeros(a.inf + sparse(I,J,complex(0,asup),m,n));
      x = intval(real(x),imag(x),'infsup');
    else
      I = find(spones(a.inf)+spones(a.sup));
      x = intval(a.inf(I),a.sup(I),'infsup');
      sizexinf = size(x.inf);
      if ( length(sizexinf)<3) & ( size(x.inf,1)==1 )
        x.inf = x.inf';
        x.sup = x.sup';
      end
    end
  end  
