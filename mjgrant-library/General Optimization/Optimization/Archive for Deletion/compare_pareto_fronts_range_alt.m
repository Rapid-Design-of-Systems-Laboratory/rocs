function compare_pareto_fronts

clear all; pack; close all; clc;


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build Pareto Fronts %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Unconstrained range vs. altitude
data = load('results_unconstrained_range_alt.mat');
in = data.in;
od = data.od;
pf_2d_unconst = build_pareto_front(od.fit(:,:,1:od.iter), ...
  'Building 05-22 Unconstrained Pareto Front');

% Constrained range vs. altitude
data = load('results_constrained_range_alt.mat');
in = data.in;
od = data.od;
pf_2d_const = build_pareto_front(od.fit(:,:,1:od.iter), ...
  'Building 05-22 Constrained Pareto Front');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unconstrained vs. Point Design %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot unconstrained Pareto front
figure;
plot(pf_2d_unconst(1,:),pf_2d_unconst(2,:),'r.');
hold on;
grid on;
% Plot 05-22 point design
plot(-7377.6,5.0538+4.1456,'k.');
legend('Unconstrained','05-22 Design');
title('05-22 Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Error Ellipse Length [km]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained vs. Unconstrained vs. Point Design %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot constrained Pareto front
figure;
plot(pf_2d_const(1,:),pf_2d_const(2,:),'b.');
hold on;
grid on;
% Plot unconstrained Pareto front
plot(pf_2d_unconst(1,:),pf_2d_unconst(2,:),'r.');
% Plot 05-22 point design
plot(-7377.6,5.0538+4.1456,'k.');
legend('Constrained','Unconstrained','05-22 Design');
axis([-7800 -6000 5.5 10]);
title('05-22 Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Error Ellipse Length [km]');

return

