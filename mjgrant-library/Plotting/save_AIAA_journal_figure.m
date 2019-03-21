function save_AIAA_journal_figure

%%%%%%%%%%%
%% Input %%
%%%%%%%%%%%

fileName = 'test';
fontSize = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alter Graph Properties %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gca,'FontSize',fontSize,'FontName','TimesNewRoman');
set(get(gca,'XLabel'),'FontSize',fontSize,'FontName','TimesNewRoman');
set(get(gca,'YLabel'),'FontSize',fontSize,'FontName','TimesNewRoman');
set(get(gca,'Title'),'FontSize',fontSize,'FontName','TimesNewRoman');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Graph with New Size and Format %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 3.25 3.25]);

% saveas(gcf,fileName,'eps');
print -depsc -r0 test.eps

return

