function vdp_plots()

clc, close all
cd data
load('results.mat')
ctrl_trad = out.setCONT(end).CONT(end).sol.control;
x_trad = out.setCONT(end).CONT(end).sol.x;

t_trad = out.setCONT(end).CONT(end).sol.parameters(1);

time_trad = t_trad*x_trad;
x1 = out.setCONT(end).CONT(end).sol.y(1,:);
x2 = out.setCONT(end).CONT(end).sol.y(2,:);
x3 = out.setCONT(end).CONT(end).sol.y(3,:);
lam2 = out.setCONT(end).CONT(end).sol.y(6,:);
lam3 = out.setCONT(end).CONT(end).sol.y(7,:);

limit = 1e-3;
ind = find(abs(lam2)<limit);
ind_sel = ind(1)

% Traditional Approach
usmooth_trad = sin(ctrl_trad(1,:));
uerror_trad = cos(ctrl_trad(1,:));
format long
lam2(end)

cd ..

for i = 1:1:length(out.setCONT(end).CONT)
    obj(i) = out.setCONT(end).CONT(i).sol.y(3,end);
    error(i) = out.setCONT(end).CONT(i).in.const.epsilon{1,1};
end

%%%%%%%%%%
%% Plot %%
%%%%%%%%%%

figure(1)
subplot(1,3,1)
yyaxis left
plot(time_trad,usmooth_trad, 'b','markersize', 3, 'linewidth', 2)
ylabel('Control, u', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
set(gca,'FontSize',12,'FontWeight' , 'bold','ycolor','b');
ylim([-1.1 1.1])
yyaxis right
plot(time_trad,lam2, 'r','markersize', 3, 'linewidth', 2);
ylabel('Switching Function, \lambda_{x_{2}}', 'fontSize', 12 , 'fontWeight' , 'bold','color','r')
ylim([-2.6 2.6])
set(gca,'FontSize',12,'FontWeight', 'bold','ycolor','r');

subplot(1,3,2)
plot(time_trad,uerror_trad, 'b','markersize', 3, 'linewidth', 2)
ylabel('Control Error, u_{\epsilon}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
set(gca,'FontSize',12,'FontWeight' , 'bold');

subplot(1,3,3)
semilogx(error,obj, 'b','markersize', 3, 'linewidth', 2)
set(gca,'FontSize',10,'FontWeight' , 'bold');
ylabel('Objective Functional, J', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('\epsilon', 'fontSize', 12 , 'fontWeight' , 'bold')

figure(2)
plot(time_trad,lam3, 'b','markersize', 3, 'linewidth', 2)
title('Switching Function History Plot', 'fontSize', 14 , 'fontWeight' , 'bold')
ylabel('\lambda_{x_{3}}', 'fontSize', 12 , 'fontWeight' , 'bold')
xlabel('Time [s]', 'fontSize', 12 , 'fontWeight' , 'bold')
grid on
set(gca,'FontSize',12,'FontWeight' , 'bold');

return