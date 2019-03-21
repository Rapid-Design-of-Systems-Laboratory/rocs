function save_fitness_func

% Open fitness_func.m
fid = fopen('fitness_func_range_alt_g.m','r');

% Read file
index = 0;
line = fgets(fid);
while line ~= -1
  index = index + 1;
  eval_obj{index} = line;
  line = fgets(fid); % Read next line
end

% Close file
fclose(fid);

data = load('results_range_alt_g.mat');
in = data.in;
od = data.od;
clear data;

% Save data
save('results_range_alt_g.mat','in','od','eval_obj');

% % Extract fitness_func.m
% data = load('test.mat');
% 
% % Check if fitness_func.m already exists
% if exist('fitness_func.m','file')
%   % File exists, error
%   error(['fitness_func.m already exists. Rename fitness_func.m to allow ', ...
%     'a new one to be created.']);
% else
%   % File does not exist, proceed with creating
%   fid = fopen('fitness_func.m','w');
% end
% 
% for ctr = 1 : 1 : length(data.eval_obj)
%   
%   fprintf(fid,'%s',data.eval_objective{ctr});
%   
% end
% 
% % Close file
% fclose(fid);

return

