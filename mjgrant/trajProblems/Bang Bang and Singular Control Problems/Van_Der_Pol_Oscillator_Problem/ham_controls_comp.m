function ham_controls_comp()
clc, close all
cd ./data
load('results_trad.mat')
ctrl = out.setCONT(end).CONT(end).sol.control;
x = out.setCONT(end).CONT(end).sol.x;
t = out.setCONT(end).CONT(end).sol.parameters(1);
time_trad = t*x;
x1_trad = out.setCONT(end).CONT(end).sol.y(1,:);
x2_trad = out.setCONT(end).CONT(end).sol.y(2,:);
lam1_trad = out.setCONT(end).CONT(end).sol.y(5,:);
lam2_trad = out.setCONT(end).CONT(end).sol.y(6,:);
lam3_trad = out.setCONT(end).CONT(end).sol.y(7,:);
lamt_trad = out.setCONT(end).CONT(end).sol.y(8,:);
usmooth_trad = sin(ctrl(1,:));
uerror_trad = cos(ctrl(1,:));
ham_trad = lam1_trad.*(x2_trad+0.001*uerror_trad)+lam2_trad.*(-x1_trad+x2_trad.*(1-x1_trad.^2)+usmooth_trad)+lam3_trad.*(0.5.*x1_trad.^2 + 0.5.*x2_trad.^2)+lamt_trad;
cd ..

cd ./data
load('results_etrig.mat')
ctrl = out.setCONT(end).CONT(end).sol.control;
x = out.setCONT(end).CONT(end).sol.x;
t = out.setCONT(end).CONT(end).sol.parameters(1);
time_trig = t*x;
x1_trig = out.setCONT(end).CONT(end).sol.y(1,:);
x2_trig = out.setCONT(end).CONT(end).sol.y(2,:);
lam1_trig = out.setCONT(end).CONT(end).sol.y(5,:);
lam2_trig = out.setCONT(end).CONT(end).sol.y(6,:);
lam3_trig = out.setCONT(end).CONT(end).sol.y(7,:);
lamt_trig = out.setCONT(end).CONT(end).sol.y(8,:);
usmooth_trig = sin(ctrl(1,:));
uerror_trig = cos(ctrl(1,:));
ham_trig = lam1_trig.*(x2_trig+0.001*uerror_trig)+lam2_trig.*(-x1_trig+x2_trig.*(1-x1_trig.^2)+usmooth_trig)+lam3_trig.*(0.5.*x1_trig.^2 + 0.5.*x2_trig.^2)+lamt_trig;
cd ..

cd ./data
load('results_pi.mat')
ctrl = out.setCONT(end).CONT(end).sol.control;
x = out.setCONT(end).CONT(end).sol.x;
t = out.setCONT(end).CONT(end).sol.parameters(1);
time_pi = t*x;
x1_pi = out.setCONT(end).CONT(end).sol.y(1,:);
x2_pi = out.setCONT(end).CONT(end).sol.y(2,:);
lam1_pi = out.setCONT(end).CONT(end).sol.y(5,:);
lam2_pi = out.setCONT(end).CONT(end).sol.y(6,:);
lam3_pi = out.setCONT(end).CONT(end).sol.y(7,:);
lamt_pi = out.setCONT(end).CONT(end).sol.y(8,:);
usmooth_pi = sin(ctrl(1,:));
uerror_pi = cos(ctrl(1,:));
ham_pi = lam1_pi.*(x2_pi+0.001*uerror_pi)+lam2_pi.*(-x1_pi+x2_pi.*(1-x1_pi.^2)+usmooth_pi)+lam3_pi.*(0.5.*x1_pi.^2 + 0.5.*x2_pi.^2)+lamt_pi;
cd ..
%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

figure(1)
subplot(1,2,1)
h1 = plot(time_trad,usmooth_trad, 'b--','markersize', 3, 'linewidth', 1.5);
ylabel('Optimal Control', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_trig,usmooth_trig, 'r*','markersize', 3, 'linewidth', 1.5);
hold on
h3 = plot(time_pi,usmooth_pi, '--','color',[0 0.5 0],'markersize', 3, 'linewidth', 1.5);
legend([h1 h2 h3],{'First Quadrant','Epsilon-Trig','Traditional'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

subplot(1,2,2)
h1 = plot(time_trad,ham_trad, 'b--','markersize', 3, 'linewidth', 1.5);
ylabel('Hamiltonian', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_trig,ham_trig, 'r*','markersize', 3, 'linewidth', 1.5);
hold on
h3 = plot(time_pi,ham_pi, '--','color',[0 0.5 0],'markersize', 3, 'linewidth', 1.5);
legend([h1 h2 h3],{'First Quadrant','Epsilon-Trig','Traditional'},'fontSize', 14)
set(gca,'FontSize',18,'FontWeight' , 'bold');

return