%-------------------------------------------------------------------------%
%                             Extract Solution                            %
%-------------------------------------------------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
state = solution.phase(1).state;
control = solution.phase(1).control;
stateExact = sin(time);

%-------------------------------------------------------------------------%
%                              Plot Solution                              %
%-------------------------------------------------------------------------%
figure(1);
pp = plot(time,state,'-o',time,stateExact,'-d');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(x(t),x^*(t))$','Interpreter','LaTeX');
ll = legend('$x(t)$','$x^*(t)$','Location','NorthWest');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'Interpreter','LaTeX','Fontsize',20);
set(gca,'Fontsize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -dpng massSpringState.png
