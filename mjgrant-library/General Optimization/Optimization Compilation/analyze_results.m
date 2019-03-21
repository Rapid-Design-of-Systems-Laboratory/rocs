function analyze_results

clear all; close all; clc;

step = 0.1;
do_plots = true;
vel_profile = [900.1 2000 5500];
bnk_profile_man = [40 40 68];

%%%%%%%%%%%%%%%%%%
%% POST Results %%
%%%%%%%%%%%%%%%%%%
% Output Best Member Info ------------------------------------------------------
bnk_profile_post = [36.8699 47.0108 50.1487];
fprintf('\n\nBest Solution - POST NPSOL');
fprintf('\n  Deployment Altitude [m]: %g',10286.3);
fprintf('\n  Early Bank Angle [deg]: %g',bnk_profile_post(1));
fprintf('\n  Mid Bank Angle [deg]: %g',bnk_profile_post(2));
fprintf('\n  Late Bank Angle [deg]: %g',bnk_profile_post(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Genetic Algorithm Results %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ga_results.mat')

% Plot contour -----------------------------------------------------------------
X = []; Y = []; fit = []; x = []; y = []; z = [];
for counter = 1 : 1 : length(ga.pop(1,1,:))
  x = [x ga.pop(1,:,counter)];
  y = [y ga.pop(2,:,counter)];
  z = [z ga.pop(3,:,counter)];
  fit = [fit ga.fitness(counter,:)];
end

if do_plots
  x_min = min(x);
  x_max = max(x);
  z_min = min(z);
  z_max = max(z);

  xg = [x_min : step : x_max];
  zg = [z_min : step : z_max];

  [X,Z] = meshgrid(xg,zg);
  FIT = griddata(x,z,fit,X,Z);

  figure;
  [c,h] = contour(X,Z,FIT,[0:1000:11000]); colorbar;
  hold on;

  % Plot actual data points
  plot(x,z,'.');

  % Plot Fitness 3-D -------------------------------------------------------------
  figure;
  plot3(x,z,fit,'.'); grid on;
  xlabel('Late Bank Angle [deg]','FontSize',14); ylabel('Early Bank Angle [deg]','FontSize',14);
  title('Deployment Altitude vs. Bank Angle Profile - GA','FontSize',14);
end

% Output Best Member Info ------------------------------------------------------
I_ga = find(fit == max(fit));
I_ga = I_ga(1); % If multiple solutions
fprintf('\n\nBest Member - Genetic Algorithm');
fprintf('\n  Deployment Altitude [m]: %g',fit(I_ga));
fprintf('\n  Early Bank Angle [deg]: %g',z(I_ga));
fprintf('\n  Mid Bank Angle [deg]: %g',y(I_ga));
fprintf('\n  Late Bank Angle [deg]: %g',x(I_ga));
bnk_profile_ga = [x(I_ga) y(I_ga) z(I_ga)];

%%%%%%%%%%%%%%%%%%%
%% Swarm Results %%
%%%%%%%%%%%%%%%%%%%
load('swarm_results.mat');

% Plot contour -----------------------------------------------------------------
X = []; Y = []; fit = []; x = []; y = []; z = [];
for counter = 1 : 1 : length(swarm.current(1,1,:))
  x = [x swarm.current(1,:,counter)];
  y = [y swarm.current(2,:,counter)];
  z = [z swarm.current(3,:,counter)];
  fit = [fit swarm.fitness(counter,:)];
end

if do_plots
  x_min = min(x);
  x_max = max(x);
  z_min = min(z);
  z_max = max(z);

  xg = [x_min : step : x_max];
  zg = [z_min : step : z_max];

  [X,Z] = meshgrid(xg,zg);
  FIT = griddata(x,z,fit,X,Z);

  figure;
  [c,h] = contour(X,Z,FIT,[0:1000:11000]); colorbar;
  hold on;

  % Plot actual data points ------------------------------------------------------
  plot(x,z,'.');

  % Plot Fitness 3-D -------------------------------------------------------------
  figure;
  plot3(x,z,fit,'.'); grid on;
  xlabel('Late Bank Angle [deg]'); ylabel('Early Bank Angle [deg]');
  title('Deployment Altitude vs. Bank Angle Profile - PSO');
end

% Output Best Member Info ------------------------------------------------------
I_swarm = find(fit == max(fit));
I_swarm = I_swarm(1); % If multiple solutions
fprintf('\n\nBest Member - PSO');
fprintf('\n  Deployment Altitude [m]: %g',fit(I_swarm));
fprintf('\n  Early Bank Angle [deg]: %g',z(I_swarm));
fprintf('\n  Mid Bank Angle [deg]: %g',y(I_swarm));
fprintf('\n  Late Bank Angle [deg]: %g\n\n',x(I_swarm));
bnk_profile_swarm = [x(I_swarm) y(I_swarm) z(I_swarm)];

% Plot profile for each best member --------------------------------------------
if do_plots
  figure; hold on; grid on;
  plot(vel_profile,bnk_profile_post,'k');
  plot(vel_profile,bnk_profile_ga,'r');
  plot(vel_profile,bnk_profile_swarm,'b');
  plot(vel_profile,bnk_profile_man,'color',[0 0.5 0]);
  legend('NPSOL','GA','PSO','Manual','Location','NorthWest');
  xlabel('Relative Velocity [m/s]','FontSize',14); ylabel('Bank Angle [deg]','FontSize',14);
%   title('Manual Optimized Bank Profile','FontSize',14);
  title('Optimized Bank Profiles','FontSize',14);
  axis([500 5500 35 75]);
end





