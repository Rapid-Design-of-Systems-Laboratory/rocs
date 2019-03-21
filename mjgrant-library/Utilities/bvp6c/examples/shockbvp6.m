function shockbvp
%SHOCKBVP  The solution has a shock layer near x = 0
%   This is an example used in U. Ascher, R. Mattheij, and R. Russell,
%   Numerical Solution of Boundary Value Problems for Ordinary Differential
%   Equations, SIAM, Philadelphia, PA, 1995,  to illustrate the mesh
%   selection strategy of COLSYS. 
%
%   For 0 < e << 1, the solution of 
%
%       e*y'' + x*y' = -e*pi^2*cos(pi*x) - pi*x*sin(pi*x)
%
%   on the interval [-1,1] with boundary conditions y(-1) = -2 and y(1) = 0
%   has a rapid transition layer at x = 0.
%
%   This example illustrates how a numerically difficult problem (e = 1e-4)
%   can be solved successfully using continuation. For this problem,
%   analytical partial derivatives are easy to derive and the solver benefits
%   from using them.  
%
%   See also BVP4C, BVPSET, BVPGET, BVPINIT, DEVAL, FUNCTION_HANDLE.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.10.4.2 $  $Date: 2005/06/21 19:27:38 $
%   BVP6C Modification
%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $

% The differential equations written as a first order system and the
% boundary conditions are coded in shockODE and shockBC, respectively. Their
% partial derivatives are coded in shockJac and shockBCJac and passed to the
% solver via the options. The option 'Vectorized' instructs the solver that
% the differential equation function has been vectorized, i.e.
% shockODE([x1 x2 ...],[y1 y2 ...]) returns [shockODE(x1,y1) shockODE(x2,y2) ...]. 
% Such coding improves the solver performance. 

% clc

tol=1e-9;
loopnum = 2;
fprintf('SHOCKBVP\n abstol = %1.1g\n reltol = %1.1g \n\n',tol,tol)

options = bvpset('reltol',tol,'abstol',tol);
options = bvpset(options,'Vectorized','on');
options = bvpset(options,'FJacobian',@shockJac,'BCJacobian',@shockBCJac);



% A guess for the initial mesh and the solution
sol = bvpinit([-1 -0.5 0 0.5 1],[1 0]);

% Problem parameter e is shared with the nested functions.
% Solution for e = 1e-2, 1e-3, 1e-4 obtained using continuation.

fprintf('  BVP4c \n')
t4 = zeros(4,1);
for l = 1:loopnum
    sol4 = sol;
    e = 0.1;
    tic
    for i=2:5
      e = e/10;   
      sol4 = bvp4c(@shockODE,@shockBC,sol4,options);
      t4(i-1) = t4(i-1)+toc;
    end
end
t4(1) = [];
disp(num2str(t4/loopnum,4))

fprintf('  BVP6c \n')
t6 = zeros(4,1);
for l = 1:loopnum
    sol6 = sol;
    e = 0.1;
    tic
    for i=2:5
      e = e/10;   
      sol6 = bvp6c(@shockODE,@shockBC,sol6,options);
      t6(i-1) = t6(i-1)+toc;
    end
end
t6(1) = [];
disp(num2str(t6/loopnum,4))

% % The final solution 
% figure;
% plot(sol.x,sol.y(1,:));
% axis([-1 1 -2.2 2.2]);
% title(['There is a shock at x = 0 when \epsilon =' sprintf('%.e',e) '.']); 
% xlabel('x');
% ylabel('solution y');

  % -----------------------------------------------------------------------
  % Nested functions -- e is shared with the outer function.
  %
  
  function dydx = shockODE(x,y)
  %SHOCKODE  Evaluate the ODE function (vectorized)
    pix = pi*x;
    dydx = [                 y(2,:)
             -x/e.*y(2,:) - pi^2*cos(pix) - pix/e.*sin(pix) ];   
  end
  % -----------------------------------------------------------------------

  function res = shockBC(ya,yb)
  %SHOCKBC  Evaluate the residual in the boundary conditions
    res = [ ya(1)+2
            yb(1)  ];
  end
  % -----------------------------------------------------------------------

  function jac = shockJac(x,y)
  %SHOCKJAC  Evaluate the Jacobian of the ODE function
  %  x and y are required arguments.
    jac = [ 0   1
            0 -x/e ];
  end
  % -----------------------------------------------------------------------

  function [dBCdya,dBCdyb] = shockBCJac(ya,yb)
  %SHOCKBCJAC  Evaluate the partial derivatives of the boundary conditions
  %  ya and yb are required arguments.
    dBCdya = [ 1 0
               0 0 ];

    dBCdyb = [ 0 0
               1 0 ];
  end
  % -----------------------------------------------------------------------

end  % shockbvp  
