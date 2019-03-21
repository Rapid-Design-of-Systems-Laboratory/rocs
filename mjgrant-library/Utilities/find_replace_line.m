function [] = find_replace_line(file,pattern,new_line,instance,offset)
%
% [] = find_replace_line(file,pattern,new_line,instance)
%
% This function is used to find and replace a line in a file specified by the 
% user.
%
% Input:
%   file - name of file [string]
%   pattern - pattern used to locate line to be modified [string]
%   new_line - new line that replaces line determined from pattern [string]
%   instance - number of times pattern must be found before replacement is made 
%               (if not assigned or empty assumes first instance) [integer]
%   offset - number of lines below pattern to perform change (useful when
%             pattern can be used to anchor which line is changed)
%
% Output:
%   None. File is modified.
%

% Author: NASA JSC DM / Michael J. Grant on March 2006


%%%%%%%%%%%%%%%
%% Open File %%
%%%%%%%%%%%%%%%

% Open file
fid = fopen(file,'r+');

% Error if file doesn't exit
if fid == -1
  error([file,' does not exist']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine instance
if ~exist('instance','var') || isempty(instance)
  instance = 1;
end

% Determine offset
if ~exist('offset','var') || isempty(offset)
  offset = 0;
end

% Read first line from file
line = fgets(fid);

% Initialize instance counter
instance_counter = 0;

% Line number - to be replaced
index_line = 1;

% Initialize index counter and text cell array
index = 1;
txt{index} = line;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Line to Replace %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read each line until match is found or end of file is reached
while line ~= -1
  
  % Check pattern
  k = strfind(line,pattern);
  
  % Determine if pattern match found
  if ~isempty(k)
    
    % Increment instance
    instance_counter = instance_counter + 1;
    
    % If correct instance, note line to replace (store in index_line)
    if instance == instance_counter
      
      % Read rest of file to cell array
      while line ~= -1
        
        % Read line
        line = fgets(fid);
        
        % Add to cell array if not end of file
        if line ~= -1
          index = index + 1;
          txt{index} = line;
        end
        
      end
      
      % End search
      break;
      
    end
    
  end
  
  % Read next line
  line = fgets(fid);
  index = index + 1;
  txt{index} = line;
  
  % Incrememt index_line
  index_line = index_line + 1;
  
end

% Close file
fclose(fid);

% Determine that a match was found
if index_line == index
  error([pattern,' not found in ',file]);
end

%%%%%%%%%%%%%%%%%%
%% Replace Line %%
%%%%%%%%%%%%%%%%%%

% Replace line in cell array and add carraige return
temp_string = txt{1};
carriage_return = temp_string(end);
txt{index_line+offset} = [new_line,carriage_return];

% Open file to be entirely replaced
fid = fopen(file,'w');

% Write new file
fprintf(fid,'%s',txt{1:end});

% Close file
fclose(fid);

return


