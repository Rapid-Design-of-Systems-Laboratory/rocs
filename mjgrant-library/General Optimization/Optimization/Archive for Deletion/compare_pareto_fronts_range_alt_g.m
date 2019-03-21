function compare_pareto_fronts

clear all; pack; close all; clc;


% Plot 3D Pareto front
figure;
plot_pareto_front('results_range_alt_g.mat','b.');
title('Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Ellipse Length [km]');
zlabel('G-Loading [Earth g]');
grid on;

% Compare 3D and 2D Range Ellipse Length vs. -Altitude
figure;
plot_pareto_front('results_range_alt_g.mat','b.',[1 2]);
hold on;
grid on;
plot_pareto_front('results_constrained_2d_range_alt_g.mat','r.');
title('Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Ellipse Length [km]');
legend('3D Pareto Surface Projection','2D Pareto Front');

% Plot 3D projection -Altitude vs. G-Loading
figure;
plot_pareto_front('results_range_alt_g.mat','b.');
title('Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Ellipse Length [km]');
zlabel('G-Loading [Earth g]');
grid on;
view(0,0);

% Plot 3D projection -Altitude vs. G-Loading
figure;
plot_pareto_front('results_range_alt_g.mat','b.');
title('Pareto Front');
xlabel('-Altitude [m]');
ylabel('Range Ellipse Length [km]');
zlabel('G-Loading [Earth g]');
grid on;
view(90,0);

return

