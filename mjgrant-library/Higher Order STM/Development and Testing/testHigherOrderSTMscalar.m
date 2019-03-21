function testHigherOrderSTMscalar

close all; clc;

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

% strEOM = 'x*sin(p(1)*x)';
% strEOM = 'sin(x)'; % This asymptotes so get interesting result
% strEOM = 'x^(1/2)'; % Shows error growing with time as expected
% strEOM = '-(x-1.1)';
% strEOM = '-x^2';
strEOM = 'x^(4/5)';
p(1) = 1;
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
tSpan = [0 2];
x0(1) = 1;
orderSTT = 4;
dx = 0.1;
tPS = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Jacobian Matrices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempJac = strEOM;
strJac = cell(orderSTT,1);
for ctr = 1 : 1 : orderSTT
  
  strJac{ctr} = diff(sym(tempJac),'x');
  tempJac = strJac{ctr};
  
  if ctr == 1
    x0(1+ctr) = 1;
  else
    x0(1+ctr) = 0;
  end
  
end

%%%%%%%%%%%%%
%% Run ODE %%
%%%%%%%%%%%%%

x0(end+1) = x0(1) + dx;
[t,x] = ode45(@eom,tSpan,x0,options,p,strEOM,strJac);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Approximate STT Solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solSTTtemp = x(:,1);
for ctr = 1 : 1 : orderSTT
  solSTT(:,ctr) = solSTTtemp + x(:,1+ctr)*dx^ctr/factorial(ctr);
  solSTTtemp = solSTT(:,ctr);
  
  errSTT(:,ctr) = solSTT(:,ctr) - x(:,end);
end

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

figure(1);
plot(t,x(:,1),'b');
hold on;
plot(t,solSTT(:,1),'r');
plot(t,solSTT(:,2),'g');
plot(t,solSTT(:,3),'c');
plot(t,solSTT(:,4),'m');
plot(t,x(:,end),'k');
grid on;
xlabel('t');
ylabel('x');
title('x vs. t');
% legend('Original Solution','Approx. STT Solution', ...
%   'Exact Perturbed Solution','Location','NorthWest');

figure(2);
plot(t,errSTT(:,1),'b');
hold on;
plot(t,errSTT(:,2),'r');
plot(t,errSTT(:,3),'c');
plot(t,errSTT(:,4),'m');
grid on;
xlabel('t');
ylabel('y error');
title('y error vs. time');

% Output partial sums at certain t
fPS = interp1(t,solSTT,tPS)

figure(3);
plot(abs(diff(fPS)));

figure(4);
plot(diff(fPS));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xdot] = eom(t,xVec,p,strEOM,strJac) %#ok<INUSL>

% Nominal solution
x = xVec(1);
xdot(1,1) = eval(strEOM);

% Compute Jacobians
for ctr = 1 : 1 : length(strJac)
  f(ctr) = eval(strJac{ctr});
end

% Compute STTs
xdot(2,1) = f(1)*xVec(2);
xdot(3,1) = f(1)*xVec(3) + f(2)*xVec(2)^2;
xdot(4,1) = f(1)*xVec(4) + f(2)*(xVec(2)*xVec(3) + xVec(3)*xVec(2) + xVec(3)*xVec(2)) + f(3)*xVec(2)^3;
xdot(5,1) = f(1)*xVec(5) + f(2)*(xVec(4)*xVec(2) + xVec(4)*xVec(2) + xVec(4)*xVec(2) + xVec(3)*xVec(3) + xVec(3)*xVec(3) + ...
  xVec(3)*xVec(3) + xVec(2)*xVec(4)) + f(3)*(xVec(3)*xVec(2)^2 + xVec(3)*xVec(2)^2 + xVec(3)*xVec(2)^2 + xVec(2)*xVec(3)*xVec(2) + ...
  xVec(2)*xVec(3)*xVec(2) + xVec(2)*xVec(2)*xVec(3)) + f(4)*xVec(2)^4;

% Exact perturbed solution
x = xVec(end);
xdot(end+1,1) = eval(strEOM);

return

