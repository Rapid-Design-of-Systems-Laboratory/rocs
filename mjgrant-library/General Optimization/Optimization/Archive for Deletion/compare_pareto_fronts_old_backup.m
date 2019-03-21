function compare_pareto_fronts_old_backup

% clear all; pack; close all; clc;

data_unconstr = load('results_unconstrained.mat'); % Unconstrained pareto front
data_constr = load('results_constrained.mat'); % Constrained Pareto front
% data_complete = load('results_constrained.mat'); % Full constrained Pareto front

figure;
% Plot unconstrained Pareto front
plot(data_unconstr.od.fit_arch(1,:,data_unconstr.od.iter), ...
  data_unconstr.od.fit_arch(2,:,data_unconstr.od.iter),'r.');
hold on;
grid on;
% Plot 05-22 point design
plot(-7377.6,5.0538+4.1456,'k.');
legend('Unconstrained','05-22 Design');
title('05-22 Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Error Ellipse Length [km]');

figure;
% Plot constrained Pareto front
plot(data_constr.od.fit_arch(1,:,data_constr.od.iter), ...
  data_constr.od.fit_arch(2,:,data_constr.od.iter),'b.');
hold on;
grid on;
% Plot unconstrained Pareto front
plot(data_unconstr.od.fit_arch(1,:,data_unconstr.od.iter), ...
  data_unconstr.od.fit_arch(2,:,data_unconstr.od.iter),'r.');
% Plot 05-22 point design
plot(-7377.6,5.0538+4.1456,'k.');
legend('Constrained','Unconstrained','05-22 Design');
axis([-7800 -6000 5.5 10]);
title('05-22 Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Error Ellipse Length [km]');

return