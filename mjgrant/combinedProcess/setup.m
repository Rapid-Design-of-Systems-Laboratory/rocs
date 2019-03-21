function setup

%%%%%%%%%%%
%% Input %%
%%%%%%%%%%%

inputDir = 'input';

%%%%%%%%%%%%%%%
%% Add Paths %%
%%%%%%%%%%%%%%%

% Add all directories of this folder to the path
new_path = genpath(pwd);
addpath(new_path);

% %%%%%%%%%%%%%%%%%%%%%%
% %% Open Input Files %%
% %%%%%%%%%%%%%%%%%%%%%%
% 
% % Get all input deck file names
% inputFileNames = dir([inputDir,'\*.m']);
% 
% % Open all input decks
% for ctr = 1 : 1 : length(inputFileNames)
%   open([inputDir,'\',inputFileNames(ctr).name]);
% end

return

