function injection
%INJECTION 
%   This example is also from U. Ascher, R. Mattheij, and R. Russell,
%   Numerical Solution of Boundary Value Problems for Ordinary Differential
%   Equations, SIAM, Philadelphia, PA, 198. here it models the injection of
%   a fluid along a long vertical channel.
%
%   The example is used here to compare performance of bvp4c and bvp6c    
%
%   See also BVP4C, BVPSET, BVPGET, BVPINIT, DEVAL, @.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.10 $  $Date: 2002/04/15 03:35:16 $
%   BVP6C Modification
%    Nick Hale  Imperial College london
%    $Date: 12/06/2006 $
% clc
tol=1e-12;
fprintf('INJECTION\n abstol = %1.1g\n reltol = %1.1g \n\n',tol,tol)

options = bvpset('stats','on','abstol',tol,'reltol',tol,'Nmax',10000);
solinit = bvpinit(linspace(0,1,10),ones(7,1),1);

fprintf('  BVP4c \n')
tic;
sol4 = bvp4c(@odes,@bcs,solinit,options);
toc

fprintf('\n  BVP6c \n')
tic;
sol6 = bvp6c(@odes,@bcs,solinit,options);
toc

figure(1)
plot(sol4.x,sol4.y(1,:),'-*',sol4.x,50.0*sol4.y(4,:),'--*', ...
    sol4.x,sol4.y(6,:),'-.*');
title(['bvp4c solutions. tol = ', num2str(tol)]);
legend('f(t)','50*h(t)','O(t)',4);
axis([0 1 0 1]);

figure(2)
plot(sol6.x,sol6.y(1,:),'-*',sol6.x,50.0*sol6.y(4,:),'--*', ...
    sol6.x,sol6.y(6,:),'-.*');
title(['bvp6c solutions. tol = ', num2str(tol)]);
legend('f(t)','50.*h(t)','O(t)',4);
axis([0 1 0 1]);

% Difference in Solutions 
x = linspace(0,1);
y4 = deval(sol4,x);
y6 = deval(sol6,x);

err=y4-y6;
for i=[1,4,6]
    Linf(i)=norm(err(i,:),inf);
end

fprintf('\n L_inf Difference in solutions of interpolation\n');
% fprintf('\ty[1] = %g\n\ty[2] = %g\n\ty[3] = %g\n\n',Linf(1),Linf(4),Linf(6));
fprintf('\t %e\n\n',max(Linf));
% -------------------------------------------------------------------------

function dydx = odes(x,y,A)
dydx = [y(2); 
        y(3);
        100*(y(2)^2-y(1)*y(3)-A);
        y(5);
        -100*y(1)*y(5)-1;
        y(7);
        -70*y(1)*y(7)];

% -------------------------------------------------------------------------

function res = bcs(ya,yb,A)
res = [ ya(1);
        ya(2);
        yb(1)-1;
        yb(2);
        ya(4);
        yb(4);
        ya(6);
        yb(6)-1 ];
