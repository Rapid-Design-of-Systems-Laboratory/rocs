function [stop options changed] = psOut(optimvalues,options,flag)

stop = false;
changed = false;
figure1 = gcf;
if strcmpi(flag,'done')
    % Create textbox
    annotation1 = annotation(...
        figure1,'textbox',...
        'Position',[0.4679 0.1357 0.3321 0.06667],...
        'FitHeightToText','off',...
        'FontWeight','bold',...
        'String',{'Pattern Search solution'});

    % Create textbox
    annotation2 = annotation(...
        figure1,'textbox',...
        'Position',[0.5625 0.2786 0.2321 0.06667],...
        'FitHeightToText','off',...
        'FontWeight','bold',...
        'String',{'FMINCON solution'});

    % Create textbox
    annotation3 = annotation(...
        figure1,'textbox',...
        'Position',[0.3714 0.6905 0.1571 0.06449],...
        'FitHeightToText','off',...
        'FontWeight','bold',...
        'String',{'Start point'});

    % Create arrow
    annotation4 = annotation(figure1,'arrow',[0.7161 0.6768],[0.3452 0.4732]);

    % Create arrow
    annotation5 = annotation(figure1,'arrow',[0.4697 0.35],[0.1673 0.2119]);

    % Create arrow
    annotation6 = annotation(figure1,'arrow',[0.4523 0.6893],[0.6929 0.6]);
end