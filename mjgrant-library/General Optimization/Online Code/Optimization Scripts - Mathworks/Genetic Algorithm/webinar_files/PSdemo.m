%% Pattern search optimization solver
% This is a demonstration of how to find a minimum of a stochastic
% objective function using PATTERNSEARCH function in the Genetic Algorithm
% and Direct Search Toolbox. We also demonstrate why The optimization
% Toolbox functions are not suitable for this kind of problems. A simple
% 2-dimensional optimization problem is selected for this demo in order to
% show some useful plots.
format compact
X0 = [2.5 -2.5];   %Starting point.
LB = [-5 -5];      %Lower bound 
UB = [5 5];        %Upper bound
range = [LB(1) UB(1); LB(2) UB(2)];
Objfcn = @smoothFcn; %Handle to the objective function.
% Plot the smooth objective function
clf;showSmoothFcn(Objfcn,range); hold on;
title('Smooth objective function')
plot3(X0(1),X0(2),feval(Objfcn,X0)+30,'om','MarkerSize',12, ...
    'MarkerFaceColor','r'); hold off;
set(gca,'CameraPosition',[-31.0391  -85.2792 -281.4265]);
set(gca,'CameraTarget',[0 0 -50])
set(gca,'CameraViewAngle',6.7937)
fig = gcf;

%% Run FMINCON on smooth objective function
% The objective function is smooth (twice continuously differentiable). We
% will solve the optimization problem using FMINCON function from the
% Optimization Toolbox. FMINCON finds a constrained minimum of a function
% of several variables. This function has a unique minimum at the point x*
% = (-5.0,-5) where it has a function value f(x*) = -250. 

% Set options to display iterative results.
options = optimset('Display','iter','OutputFcn',@fminuncOut1);
[Xop,Fop] = fmincon(Objfcn,X0,[],[],[],[],LB,UB,[],options)
figure(fig);
hold on;
%Plot the final point
plot3(Xop(1),Xop(2),Fop,'dm','MarkerSize',12,'MarkerFaceColor','m');
hold off;

%% Stochastic objective function
% The objective function is same as the previous function and some noise
% added to it.
%Reset the state of random number generators
randn('state',0);
noise = 8.5;
Objfcn = @(x) smoothFcn(x,noise); %Handle to the objective function.
%Plot the objective function (non-smooth)
figure; 
for i = 1:6
    showSmoothFcn(Objfcn,range);
    title('Stochastic objective function')
    set(gca,'CameraPosition',[-31.0391  -85.2792 -281.4265]);
    set(gca,'CameraTarget',[0 0 -50])
    set(gca,'CameraViewAngle',6.7937)
    drawnow; pause(0.2)
end
fig = gcf;
%% Run FMINCON on stochastic objective function
% The objective function is stochastic and NOT smooth. FMINCON is a general
% constrained optimization solver which finds a local minima using first
% derivative of the objective function. If derivative of the objective
% function is not provided, FMINCON uses finite difference to approximate
% first derivative  of the objective function. In this example, the
% objective function have some random noise in it. The derivatives hence
% could be highly unreliable. FMINCON can potentially stop at a point which
% is not a minimum. This may happen because the optimal conditions seems to
% be satisfied at the final point because of noise or it could not make any
% progress. 
options = optimset('Display','iter');
[Xop,Fop] = fmincon(Objfcn,X0,[],[],[],[],LB,UB,[],options)
figure(fig);
hold on;
plot3(X0(1),X0(2),feval(Objfcn,X0)+30,'om','MarkerSize',16,'MarkerFaceColor','r');
plot3(Xop(1),Xop(2),Fop,'dm','MarkerSize',12,'MarkerFaceColor','m');

%% Run PATTERNSEARCH
% We will now use PATTERNSEARCH from the Genetic Algorithm and Direct Search
% Toolbox. Pattern search optimization techniques are a class of
% direct search methods for optimization. A pattern search algorithm does
% not require any derivative information of the objective function to find
% an optimal point.
PSoptions = psoptimset('Display','iter','OutputFcn',@psOut); 
[Xps,Fps] = patternsearch(Objfcn,X0,[],[],[],[],LB,UB,PSoptions)
figure(fig);
hold on;
plot3(Xps(1),Xps(2),Fps,'pr','MarkerSize',18,'MarkerFaceColor','r');
hold off
%%
% Pattern search algorithm is not affected by random noise in the objective 
% functions. Pattern search requires only function value and not the
% derivatives, hence a noise (of some uniform kind) may not affect it.
% Pattern search requires a lot more function evaluation to find the
% minima, a cost for not using the derivatives. 

