% This is the main file for Constrained Rayleigh Problem
% input : void
% output : void
% Developed by : Kshitij Mall
% Last modified: March 18, 2016

close all; clc; clear;
nargin = 0;
if nargin < 1
    solver = 'bvp4c';
end
bvpsolver = fcnchk(solver);

Nt = 50;
tau1 = linspace(0,1,Nt); % Normalize time variable for arc 1
tau2 = 1+tau1; % Normalize time variable for arc 2
tau3 = 1+tau2; % Normalize time variable for arc 3
tau4 = 1+tau3; % Normalize time variable for arc 3
tau = [tau1 tau2 tau3 tau4]; % Normalize time variable for arc 4

% Guess for time parameters
p(1) = 1; % Time for first arc
p(2) = 0.75; % Time for second arc
p(3) = 1; % Time for third arc
p(4) = 1.75;

% Guess for unknown lambdas/nus
p(5) = 0;
p(6) = 0;
p(7) = 0;
p(8) = 0;
p(9) = 0;
p(10) = 0;

% x1 values
x1_1 = linspace(-5,-5,Nt);
x1_2 = linspace(-5,-2,Nt);
x1_3 = linspace(-2,0,Nt);
x1_4 = linspace(0,-1,Nt);

% x2 values
x2_1 = linspace(-5,5,Nt);
x2_2 = linspace(5,3,Nt);
x2_3 = linspace(3,0,Nt);
x2_4 = linspace(0,-2,Nt);

% time values
time_1 = linspace(0,1,Nt);
time_2 = linspace(0,0.75,Nt);
time_3 = linspace(0,1,Nt);
time_4 = linspace(0,1.5,Nt);

% lamx1 values
lamx1_1 = linspace(-1,1,Nt);
lamx1_2 = linspace(-1,1,Nt);
lamx1_3 = linspace(-1,1,Nt);
lamx1_4 = linspace(-1,0,Nt);

% lamx2 values
lamx2_1 = linspace(-1,-0.5,Nt);
lamx2_2 = linspace(-0.5,0.5,Nt);
lamx2_3 = linspace(0.5,0.5,Nt);
lamx2_4 = linspace(0.5,0,Nt);

% lam_time values
lam_time_1 = linspace(0.1,0.11,Nt);
lam_time_2 = linspace(0.1,0.1,Nt);
lam_time_3 = linspace(0.1,0.1,Nt);
lam_time_4 = linspace(0.1,0.1,Nt);

y1 = [x1_1;x2_1;time_1;lamx1_1;lamx2_1;lam_time_1];
y2 = [x1_2;x2_2;time_2;lamx1_2;lamx2_2;lam_time_2];
y3 = [x1_3;x2_3;time_3;lamx1_3;lamx2_3;lam_time_3];
y4 = [x1_4;x2_4;time_4;lamx1_4;lamx2_4;lam_time_4];

solinit.x = tau;
solinit.y = [y1 y2 y3 y4];
solinit.parameters = p; % Pass final time as parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
bvpOptions = bvpset('AbsTol',1e-6,'RelTol',1e-6,'Nmax',100000000,'Stats','on'); % Set options for bvp4c
tic
sol = bvpsolver(@eoms_jumps,@bcs_jumps,solinit,bvpOptions);
toc

Z = sol.y;
x1_mp = Z(1,:); % Altitude, m
x2_mp = Z(2,:); % Downrange angle, rad
time_mp = Z(3,:);
lam_x1_mp = Z(4,:); % Costate for altitude
lam_x2_mp = Z(5,:); % Costate for downrange angle
lam_t_mp = Z(6,:);
t_mp = time_mp;

for count = 1:1:length(sol.x)
    if count<126
        u_mp(count) = 1;
    elseif count>=126 && count<304
        u_mp(count) = -2*lam_x2_mp(count);
    elseif count>= 304 && count< 460
        u_mp(count) = -1;
    else
        u_mp(count) = -2*lam_x2_mp(count);
    end
end
ham_mp = lam_x1_mp.*(x2_mp)+ lam_x2_mp.*(-x1_mp + x2_mp.*(1.4 - 0.14.*x2_mp.^2) + 4.*u_mp) + lam_t_mp + u_mp.^2 + x1_mp.^2;

% Get answers for trig
load('rayleigh_trig.mat');
x1_trig = x1; % Altitude, m
x2_trig = x2; % Downrange angle, rad
time_trig = time;
lam_x1_trig = lam_x1; % Costate for altitude
lam_x2_trig = lam_x2; % Costate for downrange angle
lam_t_trig = lam_t;
t_trig = time;
u_trig = u;
ham_trig = lam_x1_trig.*x2_trig + lam_x2_trig.*(-x1_trig + x2_trig.*(1.4 - 0.14.*x2_trig.^2) + 4.*sin(u_trig)) + lam_t_trig + sin(u_trig).^2 + x1_trig.^2;

% Get answers for GPOPS-II
load('rayleigh_control.mat');
x1_gpops = x1; % Altitude, m
x2_gpops = x2; % Downrange angle, rad
time_gpops = time;
lam_x1_gpops = lamx1; % Costate for altitude
lam_x2_gpops = lamx2; % Costate for downrange angle
%lam_t_gpops = lam_t;
t_gpops = time;
u_gpops = u;
%ham_gpops = lam_x1_gpops.*x2_gpops+ lam_x2_gpops.*(-x1_gpops + x2_gpops.*(1.4 - 0.14.*x2_gpops.^2) + 4.*u_gpops) + u_gpops.^2 + x1_gpops.^2 + lam_t_gpops;
ham_gpops = ham;

ms1 = 2;
ms2 = 4;
ms3 = 2;

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

% Trajectory plot
figure(1)
subplot(2,1,1)
h1 = plot(x1_mp,x2_mp,'b', 'markersize', ms1, 'linewidth', 1.5);
xl = xlabel('x_1');
yl = ylabel('x_2');
set(xl,'FontSize',16, 'fontWeight' , 'bold');
set(yl,'FontSize',16, 'fontWeight' , 'bold');
xlim([-6.25 0.8])
set(gca,'FontSize',18,'fontWeight' , 'bold');
hold on
h2 = plot(x1_trig,x2_trig,'r-.', 'markersize', ms2, 'linewidth', 1.5);
%hold on
%h3 = plot(x1_gpops,x2_gpops,'b-+', 'markersize', ms3, 'linewidth', 1.5);
legend([h1 h2],{'Traditional','Trigonometrization'},'FontSize',14)

subplot(2,1,2)
h1 = plot(t_mp,u_mp,'b', 'markersize', ms1, 'linewidth', 1.5);
xl = xlabel('Time [s]');
yl = ylabel('u');
set(xl,'FontSize',16, 'fontWeight' , 'bold');
set(yl,'FontSize',16, 'fontWeight' , 'bold');
xlim([0 4.5])
ylim([-1.1 1.1])
set(gca,'FontSize',18,'fontWeight' , 'bold');
hold on
h2 = plot(t_trig,sin(u_trig),'r-.', 'markersize', ms2, 'linewidth', 1.5);
%h3 = plot(t_gpops,control,'b-+', 'markersize', ms3, 'linewidth', 1.5);
legend([h1 h2],{'Traditional','Trigonometrization'},'FontSize',14)

% subplot(2,2,3)
% h1 = plot(t_mp,lam_x1_mp,'r*', 'markersize', ms1, 'linewidth', 1.5);
% xl = xlabel('Time [s]');
% yl = ylabel('\lambda_{x_1}');
% set(xl,'FontSize',16, 'fontWeight' , 'bold');
% set(yl,'FontSize',16, 'fontWeight' , 'bold');
% set(gca,'FontSize',18,'fontWeight' , 'bold');
% grid on
% hold on
% h2 = plot(t_trig,lam_x1_trig,'color',[0 0.5 0],'linestyle','--', 'markersize', ms2, 'linewidth', 3);
% h3 = plot(t_gpops,lam_x1_gpops,'bx', 'markersize', ms3, 'linewidth', 1.5);
% legend([h1 h2 h3],{'Traditional','Trigonometrization', 'GPOPS-II'},'FontSize',14)
% 
% subplot(2,2,4)
% h1 = plot(t_mp,lam_x2_mp,'r*', 'markersize', ms1, 'linewidth', 1.5);
% xl = xlabel('Time [s]');
% yl = ylabel('\lambda_{x_2}');
% set(xl,'FontSize',16, 'fontWeight' , 'bold');
% set(yl,'FontSize',16, 'fontWeight' , 'bold');
% set(gca,'FontSize',18,'fontWeight' , 'bold');
% grid on
% hold on
% h2 = plot(t_trig,lam_x2_trig,'color',[0 0.5 0],'linestyle','--', 'markersize', ms2, 'linewidth', 3);
% h3 = plot(t_gpops,lam_x2_gpops,'bx', 'markersize', ms3, 'linewidth', 1.5);
% legend([h1 h2 h3],{'Traditional','Trigonometrization', 'GPOPS-II'},'FontSize',14)

figure(2)
h1 = plot(t_mp,ham_mp*1e5,'b', 'markersize', ms1, 'linewidth', 1.5);
xl = xlabel('Time [s]');
yl = ylabel('Hamiltonian X 10^{-5}');
xlim([0 4.5])
set(xl,'FontSize',24, 'fontWeight' , 'bold');
set(yl,'FontSize',24, 'fontWeight' , 'bold');
set(gca,'FontSize',20,'fontWeight' , 'bold');
hold on
h2 = plot(t_trig,ham_trig*1e5,'r-.', 'markersize', ms2, 'linewidth', 3);
%h3 = plot(t_gpops,ham_gpops,'bx', 'markersize', 3, 'linewidth', 1.5);
legend([h1 h2],{'Traditional','Trigonometrization'},'FontSize',14)

% figure(3)
% subplot(2,1,1)
% plot(t_mp,x1_mp,'r')
% title('x1 Vs time')
% xlabel('Time [s]')
% ylabel('x1')
% grid on
% hold on
% plot(t_trig,x1_trig,'b*')
% 
% subplot(2,1,2)
% plot(t_mp,x2_mp,'r')
% title('x2 Vs time')
% xlabel('Time [s]')
% ylabel('x2')
% grid on
% hold on
% plot(t_trig,x2_trig,'b*')

% End of file