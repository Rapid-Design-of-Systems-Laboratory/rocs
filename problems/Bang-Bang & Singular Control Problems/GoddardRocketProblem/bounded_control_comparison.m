function bounded_control_comparison()
close; clc
cd ./data
load('results.mat')
data1 = out.setCONT(end).CONT(end).sol;
cd ..

title_font = 20;
ylabel_font = 20;
xlabel_font = 20;
legend_font = 16;
gca_font = 20;

% UTM Result
tau_etrig = data1.x;
tf_etrig = data1.parameters(1);
time_etrig = tau_etrig*tf_etrig;
h_etrig = data1.y(1,:)/1000;
v_etrig = data1.y(2,:);
mass_etrig = data1.y(3,:);
Thrust_etrig = 96.5*sin(data1.control(1,:))+96.5;
Thrust_etrig_cos = 0.05*cos(data1.control(1,:));
Thrust_etrig_res = sqrt(Thrust_etrig.^2 + Thrust_etrig_cos.^2);

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

figure(1)

subplot(2,2,1)
h1 = plot(time_etrig,h_etrig, 'color',[0 0.5 0],'LineStyle','-','markersize', 1, 'linewidth', 5);
ylabel('Altitude [ft] X 1000', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
xlim([0 45])
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');

subplot(2,2,2)
h1 = plot(time_etrig,v_etrig, 'color',[0 0.5 0],'LineStyle','-','markersize', 1, 'linewidth', 5);
ylabel('Velocity [ft/s]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
xlim([0 45])
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');

subplot(2,2,3)
h1 = plot(time_etrig,mass_etrig, 'color',[0 0.5 0],'LineStyle','-','markersize', 1, 'linewidth', 5);
ylabel('Mass [lbm]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
xlim([0 45])
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');

subplot(2,2,4)
h1 = plot(time_etrig,Thrust_etrig, 'color',[0 0.5 0],'LineStyle','-','markersize', 1, 'linewidth', 5);
ylabel('Thrust [lbf]', 'fontSize',ylabel_font  , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
xlim([0 45])
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');

return