function make_plots

% choose t to be the tolerance you want to plot
t = 12;

tol=10^(-t);

% which results to display
one='bvp4c.solve.errors';
one='bvp5c.solve.errors';   % uncomment for bvp5c
two='bvp6c.solve.errors';

% stuff
mfc = 'markerfacecolor';
LW = 'linewidth';
lw = 1;

%%% data for bvp4c
tol4=tol;
data4=csvread([one,num2str(tol),'.txt'],2,0);
problemlist4=data4(:,1);
L2Y4=data4(:,3);
timelist4=data4(:,7);

%%% data for bvp6c
% tol = 1e-9;
tol6=tol;
data6=csvread([two,num2str(tol),'.txt'],2,0);
problemlist6=data6(:,1);
L2Y6=data6(:,3);
timelist6=data6(:,7);

% ----------- plotting --------------------

figure(1) % plot times 
plot(problemlist4,timelist4,'--k',LW,lw); hold on
plot(problemlist6,timelist6,'-k',LW,lw);
        
plot(problemlist4(L2Y4 > tol4),timelist4(L2Y4 > tol4),'sk',mfc,'w');
plot(problemlist4(L2Y4 <= tol4),timelist4(L2Y4 <= tol4),'ok',mfc,'w');
plot(problemlist6(L2Y6 > tol6),timelist6(L2Y6 > tol6),'sk',mfc,'k');
plot(problemlist6(L2Y6 <= tol6),timelist6(L2Y6 <= tol6),'ok',mfc,'k');

xlim([0,33]), set(gca,'XTick',[0:4:32])
lb = num2str([0:4:32]').'; lb(1,2:2:end) = ' '; lb(2,2:2:end) = ' '; set(gca,'XTickLabel',lb.')
xlabel('CW Problem No.','fontsize',11), ylabel('Time (secs)','fontsize',11)
text(.07,.95,['tol = 10^{-',num2str(t,2),'}'],'units','normalized','fontsize',11)



figure(2) % plot error
semilogy(problemlist4,L2Y4,'--ok',mfc,'w',LW,lw); hold on
semilogy(problemlist6,L2Y6,'-ok',mfc,'k',LW,lw); 
semilogy(problemlist6,tol,'--k');

xlim([0,33]), set(gca,'XTick',[0:4:32])
lb = num2str([0:4:32]').'; lb(1,2:2:end) = ' '; lb(2,2:2:end) = ' '; set(gca,'XTickLabel',lb.')
xlabel('CW Problem No.','fontsize',11), ylabel('L_2 error','fontsize',11)
text(.78,.95,['tol = 10^{-',num2str(t,2),'}'],'units','normalized','fontsize',11)


% % % ---------------- extra stuff for bvp4c tol = 1e-12-----------
% dataX=csvread('bvp4c_extra.txt',2,0);
% problemlistX=dataX(:,1);
% L2YX=dataX(:,3);
% timelistX=dataX(:,7);
% breakX=dataX(:,10);
% 
% figure(1), plot(problemlistX,timelistX,'^k',LW,lw);
% figure(1), semilogy(problemlistX,L2YX,'^k',mfc,'w',LW,lw); 
% % % -------------------------------------------------------------


figure(1), hold off
figure(2), hold off
