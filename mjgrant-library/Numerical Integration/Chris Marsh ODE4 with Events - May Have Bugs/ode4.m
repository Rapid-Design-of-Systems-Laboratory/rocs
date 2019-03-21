%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,x,te,xe,ie] = ode4(odefun,tspan,y0,in)
%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%
 
% if ~isnumeric(tspan)
%   error('TSPAN should be a vector of integration steps.');
% end
%  
% if ~isnumeric(y0)
%   error('Y0 should be a vector of initial conditions.');
% end
%  
h = diff(tspan);
% if any(sign(h(1))*h <= 0)
%   error('Entries of TSPAN are not in order.') 
% end  
%  
% try
%   in.prop.x0 = y0;
%   f0 = feval(odefun,tspan(1),y0,in);
% catch
%   msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
%   error(msg);  
% end  
 
y0 = y0(:);   % Make a column vector.
% if ~isequal(size(y0),size(f0))
%   error('Inconsistent sizes of Y0 and f(t0,y0).');
% end  
 
neq = length(y0);
N = length(tspan);
x = zeros(N,neq);
Y1 = zeros(neq,N);
F = zeros(neq,4);
 
Y1(:,1) = y0;
%Event checks
[value0,isterminal,direction] = Events(tspan(1),Y1(:,1),in);

terminalflag = 0;
ie = [];
te = [];
xe = [];
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  yi = Y1(:,i-1);
  in.prop.x0 = yi;
  F(:,1) = feval(odefun,ti,yi,in);
  F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),in);
  F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),in);  
  F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),in);
  Y1(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
  
  %Event checks
  [value,isterminal,direction] = Events(ti,Y1(:,i),in);
  ien = find(sign(value).*sign(value0) == -1);
  if ~isempty(ien)
      j = 1;
      while j <= length(ien)
          if (direction(ien(j)) == 1) && (sign(value0(ien(j))) > 0)
              ien(j) = [];
              j = j - 1;
          end
          if (direction(ien(j)) == -1) && (sign(value0(ien(j))) < 0)
              ien(j) = [];
              j = j - 1;
          end
          j = j + 1;
      end
  end
  
  if ~isempty(ien)
      for j = 1:length(ien)
          te = [te;ti];
          xe = [xe;Y1(:,i)'];
          if isterminal(ien(j)) == 1
              terminalflag = 1;
          end
      end
  end
  
  ie = [ie;ien];
  
  if terminalflag
      break
  end
  
  value0 = value;
  
end

t  = tspan(1:i)';
x  = Y1.';

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


