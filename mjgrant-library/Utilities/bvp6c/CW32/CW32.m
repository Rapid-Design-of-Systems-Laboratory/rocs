function CW32(problem,tol)
if nargin==0, problem=16; tol=1e-12; end
if nargin==1, tol=1e-3; end
% clc;                    close all;
fprintf('problem #%d; tol %1.1g\n',problem,tol);

%-------------------------------------SETUP--------------------------------

rtol = tol;             atol = tol;
Nstart = 33;            loopnum=1.0;

% if tol<=1e-6 loopnum=3; end
% if tol<=1e-3 loopnum=5; end

E = epslist(problem);   [xl,xr] = xlist(problem);
[ode,bc,iguess] = feval(@problemlist,problem);
if problem==31 || problem==32, dim=4; else dim=2; end; %probs 31,32 are dim 4
options = bvpset('stats','on','reltol',rtol,'abstol',atol);


%------------------------------INITIAL GUESSES-----------------------------
x=linspace(xl,xr,Nstart);

%%---zero initial guess
solinit = bvpinit(x,zeros(dim,Nstart));   
if problem==24, solinit = bvpinit(x,iguess(x,E)); end %avoids divide by zero

%%---accurate initial guess (find exact solution on perturbed mesh)
% for i=1:Nstart shift_x(i) = xl+(i-1)*(xr-xl)/(Nstart-1.0); end;
% solinit = bvpinit(x,blackbox32(problem,shift_x)); 

%%---analytic soln for linears, crude guess for non-linears
% solinit = bvpinit(x,iguess(x,E));             %


%-----------------------------------BVP4C----------------------------------

options = bvpset(options,'Nmax',floor(20000/dim));

% fprintf('\n  BVP4c \n');
tic
for i=1:loopnum
    sol4 = bvp4c(ode,bc,solinit,options,E);
end
time=toc/loopnum; 
toc

truesoln = blackbox32(problem, sol4.x);
% 
sl4x=sqrt(length(sol4.x));
Y4_L2 = norm(sol4.y(1,:)-truesoln(1,:))/sl4x;
Y4_Linf = max(abs(sol4.y(1,:)-truesoln(1,:)));
DY4_L2 = norm(sol4.y(2,:)-truesoln(2,:))/sl4x;
DY4_Linf = max(abs(sol4.y(2,:)-truesoln(2,:)));

% fprintf('bvp4c L2-Norms \n \t e(Y)  = %2.2g \n \t e(DY) = %2.2g\n',Y4_L2,DY4_L2);
% fprintf('bvp4c Linf-Norms \n \t e(Y)  = %2.2g \n \t e(DY) = %2.2g\n',Y4_Linf,DY4_Linf);

fid = fopen(['bvp4c.solve.errors' num2str(tol) '.txt'],'a+');
fprintf(fid,'%d,%0.4f,%4.4g,%4.4g,%4.4g,%4.4g,%4.4g,%d,%d,%d\n',problem,E,Y4_L2,DY4_L2,Y4_Linf,DY4_Linf,time,0,length(sol4.x),0);
fclose(fid);

%-----------------------------------BVP5C----------------------------------
if tol == 1e-12 && (problem == 32)
    options = bvpset(options,'Nmax',800);
end

ode5c = @(x,y) ode(x,y,E);
bc5c = @(x,y) bc(x,y,E);
% fprintf('\n  BVP5c \n');
tic
for i=1:loopnum
    sol5 = bvp5c(ode5c,bc5c,solinit,options);
end
time=toc/loopnum;
toc

truesoln = blackbox32(problem, sol5.x);

sl5x=sqrt(length(sol5.x));
Y5_L2 = norm(sol5.y(1,:)-truesoln(1,:))/sl5x;
Y5_Linf = max(abs(sol5.y(1,:)-truesoln(1,:)));
DY5_L2 = norm(sol5.y(2,:)-truesoln(2,:))/sl5x;
DY5_Linf = max(abs(sol5.y(2,:)-truesoln(2,:)));
% 
% fprintf('bvp5c L2-Norms \n \t e(Y)  = %2.2g \n \t e(DY) = %2.2g\n',Y5_L2,DY5_L2);
% fprintf('bvp5c Linf-Norms \n \t e(Y)  = %2.2g \n \t e(DY) = %2.2g\n',Y5_Linf,DY5_Linf);

fid = fopen(['bvp5c.solve.errors' num2str(tol) '.txt'],'a+');
fprintf(fid,'%d,%0.4f,%4.4g,%4.4g,%4.4g,%4.4g,%4.4g,%d,%d,%d\n',problem,E,Y5_L2,DY5_L2,Y5_Linf,DY5_Linf,time,0,length(sol5.x),0);
fclose(fid);


%-----------------------------------BVP6C----------------------------------

if tol == 1e-12 && (problem == 32)
    options = bvpset(options,'Nmax',200);
end

% fprintf('\n  BVP6c \n');
tic
for i=1:loopnum 
    sol6 = bvp6c(ode,bc,solinit,options,E);
end
time=toc/loopnum;
toc

truesoln = blackbox32(problem, sol6.x);

sl6x=sqrt(length(sol6.x));
Y6_L2 = norm(sol6.y(1,:)-truesoln(1,:))/sl6x;
Y6_Linf = max(abs(sol6.y(1,:)-truesoln(1,:)));
DY6_L2 = norm(sol6.y(2,:)-truesoln(2,:))/sl6x;
DY6_Linf = max(abs(sol6.y(2,:)-truesoln(2,:)));

% fprintf('bvp6c L2-Norms \n \t e(Y)  = %2.2g \n \t e(DY) = %2.2g\n',Y6_L2,DY6_L2);
% fprintf('bvp6c Linf-Norms \n \t e(Y)  = %2.2g \n \t e(DY) = %2.2g\n',Y6_Linf,DY6_Linf);

% figure(1)
% plot(sol6.x,sol6.y(1,:),'or',sol6.x,truesoln(1,:),'-b');
% figure(2)
% plot(sol6.x,sol6.y(2,:),'or',sol6.x,truesoln(2,:),'-b');

fid = fopen(['bvp6c.solve.errors' num2str(tol) '.txt'],'a+');
fprintf(fid,'%d,%0.4f,%4.4g,%4.4g,%4.4g,%4.4g,%4.4g,%d,%d,%d\n',problem,E,Y6_L2,DY6_L2,Y6_Linf,DY6_Linf,time,sol6.fevals,length(sol6.x),0);
fclose(fid);
