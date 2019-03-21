function stop = fminuncOut(x,optimvalues, state)

persistent fig gaIter
stop = false;

switch state

    case 'init'
        fig = findobj(0,'type','figure','name','Genetic Algorithm');
        limits = get(gca,'XLim');
        gaIter = limits(2);
        hold on;
    case 'iter'
        set(gca,'Xlim', [1 optimvalues.iteration + gaIter]);
        fval = optimvalues.fval;
        iter = gaIter + optimvalues.iteration;
        plot(iter,fval,'dr')
        title(['Best function value: ',num2str(fval)],'interp','none')
    case 'done'
        fval = optimvalues.fval;
        iter = gaIter + optimvalues.iteration;
        title(['Best function value: ',num2str(fval)],'interp','none')
        % Create textarrow
        annotation1 = annotation(...
            gcf,'textarrow',...
            [0.6643 0.7286],[0.3833 0.119],...
            'String',{'Algorithm switch to FMINUNC'},...
            'FontWeight','bold');
        hold off
end