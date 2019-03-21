function runCombinedProcess(inputfunc)
%
% This is the main function file that runs all the required process files
% for the optimization problem
%
% input : void
% output : void
% Developed by : Dr. M.J. Grant, Kshitij Mall and Thomas Antony
% Last modified: 25 Mar, 2014
%

% The clear all command causes segmentation fault (kernel com error) when using mathlink.
% clear all;
close all; clc;

%%%%%%%%%%%%%%%%%
%% Start Timer %%
%%%%%%%%%%%%%%%%%

tStart = tic; % Start the timer to get computational time for the entire process

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

defaultInput = defaultConfig();
% Obtain inputs for trajectory optimization
[userConfig] = inputfunc();
inTraj = mergestruct(defaultInput,userConfig);

fieldNames = fieldnames(userConfig);

for i = 1:size(fieldNames,1)
    inTraj.(fieldNames{i}) = inTraj.(fieldNames{i});
end

addpath(inTraj.autocodeDirTraj);

% Copy custom functions to autocode directory
if isfield(inTraj,'customFunc')
	for ctrFunc = 1 : 1 : length(inTraj.customFunc)
		copyfile([inTraj.customFunc(ctrFunc).name{1},'.m'],inTraj.autocodeDirTraj);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Miscellaneous Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out.origNumStates = length(inTraj.oc.state(:,1));
% out.numControl = length(inTraj.oc.control(:,1));

if inTraj.rootSolving == 1
	inTraj.convertParametersToStates = true;
end
if ~isfield(inTraj.gpuSolve,'derivativeMethod')
	inTraj.gpuSolve.derivativeMethod = 'fd';
end
% Sort constants in alphabetical order
inTraj.const = orderfields(inTraj.const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimal Control Calculations - Necessary Conditions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is used to perform all of the calculations required by indirect methods.
% The necessary conditions derived from this script are written to various
% external files to automate the process.

if true %inTraj.run.optimalCalcs

	fprintf('\nExecuting optimal control calculations...\n');

	% Unconstrained solution
	if length(inTraj.CONT) >= 1

		% Shooting does not work with parameters
		% So we always convert parameters to states when using the shooting method
		if inTraj.rootSolving == 1
			inTraj.convertParametersToStates = true;
		end

		% inTraj.oc.arcTypeSequence = 0;

		inTrajTemp = inTraj;

		% Run optimal calcs once with generic variable size inputs to account for all different problems during continuation process
		[out.oc] = optimalCalcs(inTrajTemp);
		% [out.oc, inTraj] = optimalCalcs(inTrajTemp);

		% Copy over state/costate counts
		inTraj.oc.num = out.oc.num;
		if inTraj.oc.writeEquations && inTraj.useMex
			fprintf('\nGenerating compiled code...\n');
			
      curDir = pwd;
			cd(inTraj.autocodeDirTraj);
      
      codegen computeControlUnconstrained;
      for ctrPathConstraint = 1 : 1 : inTraj.oc.num.constraints.path
        codegen(['computeControlConstraint',num2str(ctrPathConstraint)]);
      end
      
			if inTraj.rootSolving == 0
				% Create C-Mex file
				codegen derivFunc;
				codegen derivFuncRegion;
				codegen bc;
				
				if inTraj.writeJac
					codegen derivFuncJac;
				end
				% cd(curDir);
			else	% Use GPU-accelerated STM shooting
				% fprintf('Generating GPU files ...');
				% if inTraj.gpuSolve.write
					eomfile = gpuKernelSetup(@derivFuncRegion, inTraj, out);
				% else
					% eomfile = strcat(pwd,'/derivFuncRegion_stt1.cubin');
				% end
				inTraj.gpuSolve.eomfile = eomfile;
				codegen bc;

				% fprintf('done.\n');
			end
      
			cd(curDir);
		end
	end

	fprintf('done.\n');
end

% if inTraj.convertParametersToStates
% 	out.numStates = out.origNumStates + 1;
% else
% 	out.numStates = out.origNumStates;
% end
% Need modification here -- Thomas
out.numArcs   = 1;

% GPU Specific tasks
if inTraj.rootSolving == 1
	eomfile = strcat(inTraj.autocodeDirTraj,'/','derivFuncRegion','_stt',int2str(inTraj.gpuSolve.sttOrder));
	if exist(strcat(eomfile,'.cubin'),'file') > 0
		eomfile = strcat(eomfile,'.cubin');
	elseif exist(strcat(eomfile,'.ptx'),'file') > 0
		eomfile = strcat(eomfile,'.ptx');
	else
		fprintf('GPU binaries not found. Please run again with in.oc.writeEquations = true.\n');
		return;
		% eomfile = gpuKernelSetup(@derivFunc1, inTraj, inTraj.const, inTraj.constraint);
	end

	% Evaluate Jacobians of the boundary conditions
	% [inTraj.oc.M,inTraj.oc.N] = processBC(@bc,inTraj,inTraj.const, ...
	% 	[0],[],[],zeros(inTraj.oc.num.states,1),zeros(inTraj.oc.num.states,1));
	inTraj.gpuSolve.eomfile = eomfile;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Initial Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get initial guess
fprintf('\n  Constructing initial guess...');

% Make sure old format still works
if isfield(inTraj.oc,'initialGuessFunc')
    [out.IG] = inTraj.oc.initialGuessFunc(inTraj,out.oc);
else
    switch inTraj.oc.guess.mode
    case 'auto'
        [out.IG] = defaultGuess(inTraj,out.oc);
    case 'custom'
        [out.IG] = inTraj.oc.guess.func(inTraj,out.oc);
    case 'file'
        % Load from data file here
        cIndex = inTraj.oc.guess.file.index;
        iterIndex = inTraj.oc.guess.file.iteration;
        
        outBck = out;
        inBck = inTraj;
        load(inTraj.oc.guess.file.path);

        inLoad = out.setCONT(cIndex).CONT(iterIndex).in;
        sol = out.setCONT(cIndex).CONT(iterIndex).sol;
        % out = out.setCONT(cIndex).CONT(iterIndex).out;
        clear out.setCONT(cIndex+1:end);
        out.IG.sol = sol;
        inTraj.vars = inLoad.vars;
        % keyboard;

    end
end
fprintf('done.');

% <<<<<<< HEAD
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Solve Indirect Solution %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tStartBVP = tic; % Start the timer to get computational time for the indirect process
% 
% fprintf('\nExecuting Indirect Method:\n');
% if inTraj.run.indirect
% 	[out.IG] = runIndirect(inTraj,out);
% end
% 
% =======
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Execute Continuation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

tStartBVP = tic; % Start the timer to get computational time for the indirect process

fprintf('\nExecuting Continuation Method:\n');
if inTraj.run.continuation
	[inTraj, out.setCONT] = runContinuation(inTraj,out);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Final Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%

out.time = toc(tStart); % Stop timer for entire process
out.timeBVP = toc(tStartBVP); % Stop timer for continuation method
fprintf('\nContinuation process completed in %g seconds\n',out.timeBVP);
fprintf('\nFull process completed in %g seconds\n',out.time);

%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%

% Save all data
in = inTraj;
save -v7.3 data/results.mat in out;

fprintf('\nData Saved\n');

return

