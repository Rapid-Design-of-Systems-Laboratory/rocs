function bounded_control_comparison()
close all,clc
cd ./data
load('goddard_etrig.mat')
data1 = out.setCONT(end).CONT(end).sol;
load('goddard_gpopsII.mat')
time_gpopsII = time;
h_gpopsII = h;
v_gpopsII = v;
mass_gpopsII = mass;
Thrust_gpopsII = Thrust;

lamh_gpopsII = [solution.phase(1).costate(:,1);solution.phase(2).costate(:,1);solution.phase(3).costate(:,1)];
lamv_gpopsII = [solution.phase(1).costate(:,2);solution.phase(2).costate(:,2);solution.phase(3).costate(:,2)];
lammass_gpopsII = [solution.phase(1).costate(:,3);solution.phase(2).costate(:,3);solution.phase(3).costate(:,3)];
cd ..

title_font = 20;
ylabel_font = 16;
xlabel_font = 16;
legend_font = 14;
gca_font = 18;

% Epsilon Trig Result
tau_etrig = data1.x;
tf_etrig = data1.parameters(1);
time_etrig = tau_etrig*tf_etrig;
h_etrig = data1.y(1,:)/1000;
v_etrig = data1.y(2,:);
mass_etrig = data1.y(3,:);
lamh_etrig = data1.y(4,:);
lamv_etrig = data1.y(5,:);
lammass_etrig = data1.y(6,:);
Thrust_etrig = 96.5*sin(data1.control(1,:))+96.5;

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

figure(1)

subplot(2,2,1)
h1 = plot(time_etrig,h_etrig, 'bx','markersize', 4, 'linewidth', 2);
ylabel('Altitude X 1000 [ft]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpopsII,h_gpopsII, 'r--','markersize', 3, 'linewidth', 1.5);
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');
legend([h1 h2],{'Epsilon-Trig','GPOPS-II'},'fontSize',legend_font)

subplot(2,2,2)
h1 = plot(time_etrig,v_etrig/100, 'bx','markersize', 3, 'linewidth', 2);
ylabel('Velocity X 100 [ft/s]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpopsII,v_gpopsII/100, 'r--','markersize', 3, 'linewidth', 1.5);
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');
legend([h1 h2],{'Epsilon-Trig','GPOPS-II'},'fontSize',legend_font)

subplot(2,2,3)
h1 = plot(time_etrig,mass_etrig, 'bx','markersize', 3, 'linewidth', 2);
ylabel('Mass [lbm]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpopsII,mass_gpopsII, 'r--','markersize', 3, 'linewidth', 1.5);
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');
legend([h1 h2],{'Epsilon-Trig','GPOPS-II'},'fontSize',legend_font)

subplot(2,2,4)
h1 = plot(time_etrig,Thrust_etrig/100, 'bx','markersize', 3, 'linewidth', 2);
ylabel('Thrust X 100 [lbf]', 'fontSize',ylabel_font  , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpopsII,Thrust_gpopsII/100, 'r--','markersize', 3, 'linewidth', 1.5); 
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');
legend([h1 h2],{'Epsilon-Trig','GPOPS-II'},'fontSize',legend_font)

figure(2)
% Costates plot
subplot(3,1,1)
h1 = plot(time_etrig,lamh_etrig, 'bx','markersize', 4, 'linewidth', 2);
ylabel('\lambda_{h} [nd]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpopsII,lamh_gpopsII, 'r--','markersize', 3, 'linewidth', 1.5);
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');
legend([h1 h2],{'Epsilon-Trig','GPOPS-II'},'fontSize',legend_font)

subplot(3,1,2)
h1 = plot(time_etrig,lamv_etrig/10, 'bx','markersize', 3, 'linewidth', 2);
ylabel('\lambda_{v} X 10 [s]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpopsII,lamv_gpopsII/10, 'r--','markersize', 3, 'linewidth', 1.5);
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');
legend([h1 h2],{'Epsilon-Trig','GPOPS-II'},'fontSize',legend_font)

subplot(3,1,3)
h1 = plot(time_etrig,lammass_etrig/10000, 'bx','markersize', 3, 'linewidth', 2);
ylabel('\lambda_{mass} X 10^{4} [ft/lbm]', 'fontSize', ylabel_font , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', xlabel_font , 'fontWeight' , 'bold')
grid on
hold on
h2 = plot(time_gpopsII,lammass_gpopsII/10000, 'r--','markersize', 3, 'linewidth', 1.5);
set(gca,'FontSize',gca_font,'FontWeight' , 'bold');
legend([h1 h2],{'Epsilon-Trig','GPOPS-II'},'fontSize',legend_font)
return