%----------------------------------%
%   Plot Solution from GPOPS-II    %
%----------------------------------%
solution = output.result.solution;
marks = {'-bd','-gs','-r^','-cv','-mo','-kh'};
% Extract the time in each phase of the problem.
time{1} = solution.phase(1).time;
time{2} = solution.phase(2).time;
time{3} = solution.phase(3).time;
timeRadau{1} = solution.phase(1).timeRadau;
timeRadau{2} = solution.phase(2).timeRadau;
timeRadau{3} = solution.phase(3).timeRadau;
for iphase=1:3
  y{iphase} = solution.phase(iphase).state;
  u{iphase} = solution.phase(iphase).control;
  uRadau{iphase} = solution.phase(iphase).controlRadau;
end;
nstates = 6;
ncontrols = 5;
for istate = 1:nstates
  stateString = strcat('$y_',num2str(istate),'(t)$');
  stateStringPrint = strcat('Y',num2str(istate));
  figure(istate);
  pp = plot(time{1},y{1}(:,istate),marks{1},time{2},y{2}(:,istate),marks{2},time{3},y{3}(:,istate),marks{3});
  xl = xlabel('$t \textrm{(hours)}$','interpreter','LaTeX');
  yl = ylabel(stateString,'interpreter','LaTeX');
  ll = legend('$\textrm{Phase 1}$','$\textrm{Phase 2}$','$\textrm{Phase 3}$','Location','Best');
  set(pp,'LineWidth',1.25,'MarkerSize',7);
  set(xl,'FontSize',18,'interpreter','LaTeX');
  set(yl,'FontSize',18,'interpreter','LaTeX');
  set(ll,'FontSize',18,'interpreter','LaTeX');
  set(gca,'FontSize',16,'FontName','Times');
  grid on;
  filename = strcat([,'kineticBatchReactor',stateStringPrint,'.eps']);
  print('-depsc2',filename);
  filename = strcat([,'kineticBatchReactor',stateStringPrint,'.png']);
  print('-dpng',filename);
end;

for icontrol = 1:ncontrols
  figure(icontrol+nstates);
  controlStringPrint = strcat('U',num2str(icontrol));
  if (icontrol==2),
    controlString = strcat('$u_',num2str(icontrol),'(t)\times 10^{-3}$');
    pp = plot(timeRadau{1},uRadau{1}(:,icontrol)*1e3,marks{1},timeRadau{2},uRadau{2}(:,icontrol)*1e3,marks{2},timeRadau{3},uRadau{3}(:,icontrol)*1e3,marks{3});  
  elseif (icontrol==3) || (icontrol==4),
    controlString = strcat('$u_',num2str(icontrol),'(t)\times 10^{-5}$');
    pp = plot(timeRadau{1},uRadau{1}(:,icontrol)*1e5,marks{1},timeRadau{2},uRadau{2}(:,icontrol)*1e5,marks{2},timeRadau{3},uRadau{3}(:,icontrol)*1e5,marks{3});
  else
    controlString = strcat('$u_',num2str(icontrol),'(t)$');
    pp = plot(timeRadau{1},uRadau{1}(:,icontrol),marks{1},timeRadau{2},uRadau{2}(:,icontrol),marks{2},timeRadau{3},uRadau{3}(:,icontrol),marks{3});
  end
  xl = xlabel('$t \textrm{(hr)}$','interpreter','LaTeX');
  yl = ylabel(controlString,'interpreter','LaTeX');
  ll = legend('$\textrm{Phase 1}$','$\textrm{Phase 2}$','$\textrm{Phase 3}$','Location','Best');
  set(pp,'LineWidth',1.25,'MarkerSize',7);
  set(xl,'FontSize',18,'interpreter','LaTeX');
  set(yl,'FontSize',18,'interpreter','LaTeX');
  set(ll,'FontSize',18,'interpreter','LaTeX');
  set(gca,'FontSize',16,'FontName','Times');
  grid on;
  filename = strcat([,'kineticBatchReactor',controlStringPrint,'.eps']);
  print('-depsc2',filename);
  filename = strcat([,'kineticBatchReactor',controlStringPrint,'.png']);
  print('-dpng',filename);
end;

n = nstates + ncontrols;

for istate = 1:nstates
  stateString = strcat('$y_',num2str(istate),'(t)$');
  stateStringPrint = strcat('Y',num2str(istate));
  figure(istate + n);
  pp = plot(time{1},y{1}(:,istate),marks{1});
  xl = xlabel('$t \textrm{(hours)}$','interpreter','LaTeX');
  yl = ylabel(stateString,'interpreter','LaTeX');
  ll = legend('$\textrm{Phase 1}$');
  set(pp,'LineWidth',1.25,'MarkerSize',7);
  set(xl,'FontSize',18,'interpreter','LaTeX');
  set(yl,'FontSize',18,'interpreter','LaTeX');
  set(ll,'FontSize',18,'interpreter','LaTeX');
  set(gca,'FontSize',16,'FontName','Times');
  grid on;
  %filename = strcat([,'kineticBatchReactor',stateStringPrint,'Phase1.eps']);
  %print('-depsc2',filename);
end;

for icontrol = 1:ncontrols
  figure(icontrol+nstates + n);
  controlStringPrint = strcat('U',num2str(icontrol));
  if (icontrol==2),
    controlString = strcat('$u_',num2str(icontrol),'(t)\times 10^{-3}$');
    pp = plot(timeRadau{1},uRadau{1}(:,icontrol)*1e3,marks{1});
  elseif (icontrol==3) || (icontrol==4),
    controlString = strcat('$u_',num2str(icontrol),'(t)\times 10^{-5}$');
    pp = plot(timeRadau{1},uRadau{1}(:,icontrol)*1e5,marks{1});
  else
    controlString = strcat('$u_',num2str(icontrol),'(t)$');
    pp = plot(timeRadau{1},uRadau{1}(:,icontrol),marks{1});
  end
  xl = xlabel('$t \textrm{(hours)}$','interpreter','LaTeX');
  yl = ylabel(controlString,'interpreter','LaTeX');
  ll = legend('$\textrm{Phase 1}$');
  set(pp,'LineWidth',1.25,'MarkerSize',7);
  set(xl,'FontSize',18,'interpreter','LaTeX');
  set(yl,'FontSize',18,'interpreter','LaTeX');
  set(ll,'FontSize',18,'interpreter','LaTeX');
  set(gca,'FontSize',16,'FontName','Times');
  grid on;
  %filename = strcat([,'kineticBatchReactor',controlStringPrint,'Phase1.eps']);
  %print('-depsc2',filename);
end;
