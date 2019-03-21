function [] = compile_optim()
%
% [] = compile_optim()
%
% This function is used to compile the optimization routines into an
% executable. This executable does not require a Matlab license to run.
% The executable is called ea_optim. A bash scrip (run_optim) is also
% created. run_optim sets all of the necessary settings by pointing Linux
% to the Matlab headerless library installed in the home directory.
% run_optim then executes ea_optim. Thus, the optimization is executed
% by simply typing 'run_optim' into the command line (without the quotes).
% Note that the optimization does not seem to work if run_optim is
% executed in the background.
%
% Input:
%   None.
%
% Output:
%   None. However, several files are created with the executable.
%

% Author: NASA JSC/DM42 - Michael J. Grant


%%%%%%%%%%%%%
%% Compile %%
%%%%%%%%%%%%%

% Set path information. Assumes path to optimization directory added by
% set_input.m
[in] = set_input;

% Compile optimization algorithm
eval('mcc -m ea_optim');

% Create executable script
fid = fopen('run_optim','w');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'\n');
fprintf(fid,'ROOT=/misc/home0/mjgrant/matlab/MCR/v73\n');
fprintf(fid,'\n');
fprintf(fid,'export MATLAB=$ROOT\n');
fprintf(fid,'\n');
fprintf(fid,'export LD_LIBRARY_PATH=$ROOT/runtime/glnx86:\\\n');
fprintf(fid,'$ROOT/bin/glnx86:\\\n');
fprintf(fid,'$ROOT/sys/os/glnx86:\\\n');
fprintf(fid,'$ROOT/sys/java/jre/glnx86/jre1.5.0/lib/i386/client:\\\n');
fprintf(fid,'$ROOT/sys/java/jre/glnx86/jre1.5.0/lib/i386:\\\n');
fprintf(fid,'$ROOT/sys/opengl/lib/glnx86:LD_LIBRARY_PATH\n');
fprintf(fid,'\n');
fprintf(fid,'export XAPPLRESDIR=$ROOT/X11/app-defaults\n');
fprintf(fid,'\n');
fprintf(fid,'./ea_optim\n');
fprintf(fid,'\n');
fclose(fid);

!chmod u+x run_optim

return


