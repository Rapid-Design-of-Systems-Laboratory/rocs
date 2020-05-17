function comparison_brach()
clc, close all
cd ./data

% Traditional Way
load('modbrachapprox.mat')
data1 = out.setCONT(end).CONT(end).sol;
alfa_approx = data1.control;
time_approx = data1.x;
tf_approx = data1.parameters(1);
t_approx = tf_approx*time_approx;
x_approx = data1.y(1,:);
y_approx = data1.y(2,:);
v_approx = data1.y(3,:);

for count = 2:1:length(alfa_approx)
    time_fac = t_approx(1,1:count);
    alfa_approx_fac = alfa_approx(1,1:count)*180/pi;
    cost_approx(1,count) = trapz(time_fac,alfa_approx_fac.^2/2);
end

% Trig Way
load('modbrachtrig.mat')
data2 = out.setCONT(end).CONT(end).sol;
time_trig = data2.x;
tf_trig = data2.parameters(1);
t_trig = tf_trig*time_trig;
x_trig = data2.y(1,:);
y_trig = data2.y(2,:);
v_trig = data2.y(3,:);
alfa_trig = data2.y(4,:);

for count = 2:1:length(alfa_trig)
    time_fac = t_trig(1,1:count);
    alfa_trig_fac = alfa_trig(1,1:count)*180/pi;
    cost_trig(1,count) = trapz(time_fac,alfa_trig_fac.^2/2);
end

% GPOPS-II 
load('modbrachgpops.mat')
data3 = solution.phase(end);
alfa_gpops2 = data3.control;
t_gpops2 = data3.time;
x_gpops2 = data3.state(:,1);
y_gpops2 = data3.state(:,2);
v_gpops2 = data3.state(:,3);

for count = 2:1:length(alfa_gpops2)
    time_fac = t_gpops2(1:count);
    alfa_gpops2_fac = alfa_gpops2(1:count)*180/pi;
    cost_gpops2(1,count) = trapz(time_fac,alfa_gpops2_fac.^2/2);
end

cd ..

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

figure(1)
subplot(2,2,1)
h1 = plot(t_approx,x_approx, 'r--','markersize', 2, 'linewidth', 2);
ylabel('x [m]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
hold on
h2 = plot(t_trig,x_trig, 'b*','markersize', 2, 'linewidth', 2);
h3 = plot(t_gpops2,x_gpops2, 'color',[0 0.5 0],'LineStyle','--','markersize', 1, 'linewidth', 1);
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2 h3],{'Approximation','Control-Rate','GPOPS-II',},'fontSize',12)
xlim([0 0.83])

subplot(2,2,2)
h1 = plot(t_approx,y_approx, 'r--','markersize', 2, 'linewidth', 2);
ylabel('y [m]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16  , 'fontWeight' , 'bold')
hold on
h2 = plot(t_trig,y_trig, 'b*','markersize', 2, 'linewidth', 2);
h3 = plot(t_gpops2,y_gpops2, 'color',[0 0.5 0],'LineStyle','--','markersize', 1, 'linewidth', 1);
set(gca,'FontSize',16 ,'FontWeight' , 'bold');
legend([h1 h2 h3],{'Approximation','Control-Rate','GPOPS-II',},'fontSize',12)
xlim([0 0.83])

subplot(2,2,3)
h1 = plot(t_approx,v_approx, 'r--','markersize', 2, 'linewidth', 2);
ylabel('v [m/s]', 'fontSize', 16  , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16  , 'fontWeight' , 'bold')
hold on
h2 = plot(t_trig,v_trig, 'b*','markersize', 2, 'linewidth', 2);
h3 = plot(t_gpops2,v_gpops2, 'color',[0 0.5 0],'LineStyle','--','markersize', 1, 'linewidth', 1);
set(gca,'FontSize',16 ,'FontWeight' , 'bold');
legend([h1 h2 h3],{'Approximation','Control-Rate','GPOPS-II',},'fontSize',12)
xlim([0 0.83])

subplot(2,2,4)
h1 = plot(t_approx,alfa_approx*180/pi, 'r--','markersize', 2, 'linewidth', 2);
ylabel('\alpha [deg]', 'fontSize', 16  , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16  , 'fontWeight' , 'bold')
hold on
h2 = plot(t_trig,alfa_trig*180/pi, 'b*','markersize', 2, 'linewidth', 2);
h3 = plot(t_gpops2,alfa_gpops2*180/pi, 'color',[0 0.5 0],'markersize', 1, 'linewidth', 1);
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2 h3],{'Approximation','Control-Rate','GPOPS-II',},'fontSize',12)
xlim([0 0.83])
ylim([0 81])

figure(2)
h1 = plot(t_approx,cost_approx, 'r--','markersize', 2, 'linewidth', 2);
ylabel('Cost Functional [deg^2/s]', 'fontSize', 16 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 16 , 'fontWeight' , 'bold')
hold on
h2 = plot(t_trig,cost_trig, 'b*','markersize', 2, 'linewidth', 2);
h3 = plot(t_gpops2,cost_gpops2, 'color',[0 0.5 0],'LineStyle','--','markersize', 1, 'linewidth', 1);
set(gca,'FontSize',16,'FontWeight' , 'bold');
legend([h1 h2 h3],{'Approximation','Control-Rate','GPOPS-II',},'fontSize',12)
xlim([0 0.83])
return