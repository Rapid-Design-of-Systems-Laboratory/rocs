function create_folder(folder_name)
%
% This function is used to delete an existing folder (and all of its contents),
% if it exists, and re-create the folder (so folder is empty).
%

% Author: GaTech SSDL - Michael J. Grant on 11/01/08

% Check if folder exists. If it does, delete it.
if exist(folder_name,'dir')
  rmdir(folder_name,'s');
end

% Create new folder
mkdir(folder_name);

return

