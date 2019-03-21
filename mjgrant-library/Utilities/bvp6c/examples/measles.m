function measles
%MEASLES 
%   This example is from U. Ascher, R. Mattheij, and R. Russell,
%   Numerical Solution of Boundary Value Problems for Ordinary Differential
%   Equations, SIAM, Philadelphia, PA, 1988 which models the spread of
%   measles in a population.
%
%   The example is used here to compare performance of bvp4c and bvp6c    
%
%   See also BVP4C, BVPSET, BVPGET, BVPINIT, DEVAL, @.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.10 $  $Date: 2002/04/15 03:35:16 $
%   BVP6C Modification
%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $
tol=1e-6;
fprintf('MEASLES\n abstol = %1.1g\n reltol = %1.1g \n\n',tol,tol)

options = bvpset('Stats','on','reltol',tol,'abstol',tol);
solinit = bvpinit(linspace(0,1,16),[0.01, 0.01, 0.01]);

fprintf('  BVP4c \n')
tic;
sol4 = bvp4c(@odes,@bcs,solinit,options);
toc

options = bvpset('ls
Stats','on','reltol',tol,'abstol',tol);
fprintf('\n  BVP6c \n')
tic;
sol6 = bvp6c(@odes,@bcs,solinit,options);
toc

% figure(1)
% plot(sol4.x,sol4.y(1,:)-0.07,'-',sol4.x,sol4.y(2,:),'--', ...
%     sol4.x,sol4.y(3,:),'-.');
% title(['bvp4c solutions. tol = ', num2str(tol)]);
% legend('y_1(t) - 0.07 ','y_2(t)','y_3(t)',4);
% 
% figure(2)
% plot(sol6.x,sol6.y(1,:)-0.07,'-',sol6.x,sol6.y(2,:),'--', ...
%     sol6.x,sol6.y(3,:),'-.');
% title(['bvp6c solutions. tol = ', num2str(tol)]);
% legend('y_1(t) - 0.07 ','y_2(t)','y_3(t)',4);

% Difference in Solutions 
x = linspace(0,1);

y4 = deval(sol4,x);
y6 = deval(sol6,x);

err=y4-y6;
for i=1:3
    Linf(i)=norm(err(i,:),inf);
end

fprintf('\n L_inf Difference in solutions of interpolation\n');
% fprintf('\ty[1] = %g\n\ty[2] = %g\n\ty[3] = %g\n\n',Linf(1),Linf(2),Linf(3));
fprintf('\t %e\n\n',max(Linf));

% -------------------------------------------------------------------------
function dydx = odes(x,y)
beta = 1575*(1+cos(2*pi*x));
dydx = [    0.02 - beta*y(1)*y(3)
        beta*y(1)*y(3) - y(2)/0.0279
        y(2)/0.0279 - y(3)/0.01     ];
    


% -------------------------------------------------------------------------

function res = bcs(ya,yb)
res = ya - yb;
