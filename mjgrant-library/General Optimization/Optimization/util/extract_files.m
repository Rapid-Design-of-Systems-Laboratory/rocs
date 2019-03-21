function [] = extract_files(file)
%
% [] = extract_files(file)
%
% This function is used to extract fitness_func.m and set_input.m from the
% results file created during the optimization.
%
% Input:
%   file - file name of optimization results file [string]
%
% Output:
%   None (fitness_func.m and set_input.m are extracted)
%

% Author: NASA JSC/DM42 - Michael J. Grant


%%%%%%%%%%%%%%%%%
%% Input Logic %%
%%%%%%%%%%%%%%%%%

if ~exist('file','var')
  file = 'results.mat';
end

% Extract fitness_func.m
temp = load(file,'ui');
ui = temp.ui;
clear temp;

%%%%%%%%%%%%%%%%%
%% Error Check %%
%%%%%%%%%%%%%%%%%

% Check if fitness_func.m already exists
if exist('./eval_pop.m','file')
  % File exists, error
  error(['./fitness_func.m already exists. Rename fitness_func.m to allow ', ...
    'a new one to be created.']);
end

% Check if set_input.m already exists
if exist('./set_input.m','file')
  % File exists, error
  error(['./set_input.m already exists. Rename set_input.m to allow ', ...
    'a new one to be created.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract fitness_func.m %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open fitness_func.m for writing
fid = fopen('./eval_pop.m','w');

% Write to fitness_func.m
for ctr = 1 : 1 : length(ui.eval_obj)
  fprintf(fid,'%s',ui.eval_obj{ctr});
end

% Close file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract set_input.m %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Open fitness_func.m for writing
fid = fopen('./set_input.m','w');

% Write to fitness_func.m
for ctr = 1 : 1 : length(ui.set_input)
  fprintf(fid,'%s',ui.set_input{ctr});
end

% Close file
fclose(fid);

return

