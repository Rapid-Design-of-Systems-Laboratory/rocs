function testAccuracyTS

close all; clc;

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

f{1} = 'log(1+x)'; % Function
plotSet = [3 6 9 15];
% f{1} = 'exp(x)'; % Function - good example of convergence without oscillation
% plotSet = [2 3 4 6];
% f{1} = 'cos(x)+sin(x)'; % Function
a = 0; % Reference point for Taylor Series approximation
n = 21; % Number of Taylor Series terms
x = linspace(-1,1.4,61)';
plotColor = {'k:','k-.','k--','k'};
xPS = [0.5 0.9 1.1 1.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Taylor Series Terms and Construct Solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fVal{1} = subs(sym(f{1}),sym('x'),a);
solApprox{1} = sym(fVal{1});

for ctr = 1 : 1 : n-1
  
  % Compute Jacobian term
  f{ctr+1} = diff(sym(f{ctr}),'x');
  fVal{ctr+1} = subs(sym(f{ctr+1}),sym('x'),a);
  
  % Compute approximate solution
  solApprox{ctr+1} = solApprox{ctr} + ...
    sym(fVal{ctr+1})/factorial(ctr)*(sym('x')-a)^(ctr);
  
end

solExact = eval(sym(f{1}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Taylor Series Solutions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ctr = 1 : 1 : n
  
  % Check if answer independent of x. If so, need to repeat value to make
  % same size in matrix as other terms.
  tempVal = eval(solApprox{ctr});
  if length(tempVal) == 1
    tempVal = tempVal*ones(length(x),1);
  end
  
  solApproxVal(:,ctr) = tempVal;
  
end

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

% Plot solutions
figure(1); hold on; grid on;
xlabel('x');
ylabel(f{1});

index = 0;
for ctr = plotSet
  
  index = index + 1;
  plot(x,solApproxVal(:,ctr),plotColor{index},'LineWidth',1);
  
end

plot(x(1:3:end),solExact(1:3:end),'k*','LineWidth',1);
legendStr = ['legend(',sprintf('''p = %i'',',plotSet),'''Exact'',''Location'',''NorthWest'');'];
eval(legendStr);
presentation_plot;

% Output partial sums at certain x
for ctr = 1 : 1 : length(xPS)
  fPS{ctr} = interp1(x,solApproxVal,xPS(ctr));
end

for ctr = 1 : 1 : length(xPS)
  
  figure(2);
  semilogy(abs(diff(fPS{ctr})),plotColor{ctr},'LineWidth',1);
%   plot(abs(diff(fPS{ctr})),plotColor{ctr},'LineWidth',1);
  xlabel('Order p');
  strTemp = f{1};
  ylabel('|f_p_+_1 - f_p|');
  hold on;
  grid on;

  figure(3);
  plot(diff(fPS{ctr}),plotColor{ctr},'LineWidth',1);
  hold on;
  presentation_plot;

end

figure(2);
legendStr = ['legend(',sprintf('''x = %g'',',xPS),'''Location'',''NorthEast'');'];
eval(legendStr);
presentation_plot;

return

