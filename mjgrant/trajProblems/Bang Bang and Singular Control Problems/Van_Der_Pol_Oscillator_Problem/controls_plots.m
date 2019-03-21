function controls_new()
clc, close all
cd ./data
load('results_etrig.mat')
% load('results.mat')
ctrl_trad = out.setCONT(end).CONT(end).sol.control;
x_trad = out.setCONT(end).CONT(end).sol.x;

t_trad = out.setCONT(end).CONT(end).sol.parameters(1);

time_trad = t_trad*x_trad;
x1 = out.setCONT(end).CONT(end).sol.y(1,:);
x2 = out.setCONT(end).CONT(end).sol.y(2,:);
x3 = out.setCONT(end).CONT(end).sol.y(3,:);
lam1 = out.setCONT(end).CONT(end).sol.y(5,:);
lam2 = out.setCONT(end).CONT(end).sol.y(6,:);
lam3 = out.setCONT(end).CONT(end).sol.y(7,:);
lamt = out.setCONT(end).CONT(end).sol.y(8,:);
cd ..
% limit = 0.003;
% ind = find(abs(lam2)<limit);
ind_sel = 135;

% Traditional Approach
usmooth_trad = sin(ctrl_trad(1,:));
uerror_trad = cos(ctrl_trad(1,:));

% Non singular solution
time_non = time_trad(1:ind_sel);
x1_non = x1(1:ind_sel);
x2_non = x2(1:ind_sel);
x3_non = x3(1:ind_sel);
lam1_non = lam1(1:ind_sel);
lam2_non = lam2(1:ind_sel);
lam3_non = lam3(1:ind_sel);
lamt_non = lamt(1:ind_sel);
usmooth_non = usmooth_trad(1:ind_sel);
uerror_non = uerror_trad(1:ind_sel);

% Near singular solution
time_near = time_trad(ind_sel:end);
x1_near = x1(ind_sel:end);
x2_near = x2(ind_sel:end);
x3_near = x3(ind_sel:end);
lam1_near = lam1(ind_sel:end);
lam2_near = lam2(ind_sel:end);
lam3_near = lam3(ind_sel:end);
lamt_near = lamt(ind_sel:end);
usmooth_near = usmooth_trad(ind_sel:end);
uerror_near = uerror_trad(ind_sel:end);

ham_non = lam1_non.*(x2_non+0.001*uerror_non)+lam2_non.*(-x1_non+x2_non.*(1-x1_non.^2)+usmooth_non)+lam3_non.*(0.5.*x1_non.^2 + 0.5.*x2_non.^2)+lamt_non;
ham_near = lam1_near.*(x2_near+0.001*uerror_near)+lam2_near.*(-x1_near+x2_near.*(1-x1_near.^2)+usmooth_near)+lam3_near.*(0.5.*x1_near.^2 + 0.5.*x2_near.^2)+lamt_near;
%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

figure(1)
subplot(1,2,1)
h1 = plot(time_non,usmooth_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('u', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near,usmooth_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(1,2,2)
h1 = plot(time_non,uerror_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('u_{\epsilon}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near,uerror_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');


figure(2)
subplot(3,2,1)
h1 = plot(time_non, x1_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('x_{1}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near, x1_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(3,2,2)
h1 = plot(time_non, x2_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('x_{2}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near, x2_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(3,2,3)
h1 = plot(time_non, x3_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('x_{3}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near, x3_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(3,2,4)
h1 = plot(time_non,lam1_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('\lambda_{x_{1}}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near,lam1_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(3,2,5)
h1 = plot(time_non,lam2_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('\lambda_{x_{2}}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near,lam2_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(3,2,6)
h1 = plot(time_non,lam3_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('\lambda_{x_{3}}', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near,lam3_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

figure(3)
h1 = plot(time_non,ham_non, 'b','markersize', 3, 'linewidth', 2);
ylabel('ham', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_near,ham_near, 'r','markersize', 3, 'linewidth', 2);
legend([h1 h2],{'Bang-Bang','Near-Singular'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');
return