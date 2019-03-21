figure(1)
pp = plot(solution.phase(1).time,solution.phase(1).state,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(x(t),y(t),v(t))$','Interpreter','LaTeX');
ll = legend('$x(t)$','$y(t)$','$v(t)$','Location','NorthWest');
set(pp,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(gca,'FontSize',16,'FontName','Times');
grid on
print -dpng brachistochroneState.png

figure(2)
pp = plot(solution.phase(1).time,solution.phase(1).control,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
ll = legend('$u(t)$','Location','NorthWest');
set(pp,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(gca,'FontSize',16,'FontName','Times');
grid on
print -dpng brachistochroneControl.png
