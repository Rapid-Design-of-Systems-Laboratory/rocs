% gpops2PathSetup
% this script adds the appropriate paths for use of gpops2 to the
% MATLAB path directory

% get current directory
currdir = pwd;

% add adigator directories
% if isdir('adigator'),
%   disp(['Adding Directory ',currdir,'/adigator/lib/ to Path']);
%   addpath([currdir,'/adigator/lib/'],'-begin');
%   
%   disp(['Adding Directory ',currdir,'/adigator/util/ to Path']);
%   addpath([currdir,'/adigator/util/'],'-begin');
% end;

% add RPMintegration/opRPMIsnopt/
disp(['Adding Directory ',currdir,'/lib/gpopsRPMIntegration/gpopsSnoptRPMI/ to Path']);
addpath([currdir,'/lib/gpopsRPMIntegration/gpopsSnoptRPMI/'],'-begin');

% add RPMintegration/opRPMIipopt/
disp(['Adding Directory ',currdir,'/lib/gpopsRPMIntegration/gpopsIpoptRPMI/ to Path']);
addpath([currdir,'/lib/gpopsRPMIntegration/gpopsIpoptRPMI/'],'-begin');

% add RPMintegration/
disp(['Adding Directory ',currdir,'/lib/gpopsRPMIntegration/ to Path']);
addpath([currdir,'/lib/gpopsRPMIntegration/'],'-begin');

% add RPMdifferentiation/opRPMDsnopt/
disp(['Adding Directory ',currdir,'/lib/gpopsRPMDifferentiation/gpopsSnoptRPMD/ to Path']);
addpath([currdir,'/lib/gpopsRPMDifferentiation/gpopsSnoptRPMD/'],'-begin');

% add RPMdifferentiation/opRPMDipopt/
disp(['Adding Directory ',currdir,'/lib/gpopsRPMDifferentiation/gpopsIpoptRPMD/ to Path']);
addpath([currdir,'/lib/gpopsRPMDifferentiation/gpopsIpoptRPMD/'],'-begin');

% add RPMdifferentiation/
disp(['Adding Directory ',currdir,'/lib/gpopsRPMDifferentiation/ to Path']);
addpath([currdir,'/lib/gpopsRPMDifferentiation/'],'-begin');

% add OCPfinitediff/
disp(['Adding Directory ',currdir,'/lib/gpopsFiniteDifference/ to Path']);
addpath([currdir,'/lib/gpopsFiniteDifference/'],'-begin');

% add AdiGator/
disp(['Adding Directory ',currdir,'/lib/gpopsADiGator/ to Path']);
addpath([currdir,'/lib/gpopsADiGator/'],'-begin');

% add gpopsAutomaticScaling/
disp(['Adding Directory ',currdir,'/lib/gpopsAutomaticScaling/ to Path']);
addpath([currdir,'/lib/gpopsAutomaticScaling/'],'-begin');

% add gpopsMeshRefinement/
disp(['Adding Directory ',currdir,'/lib/gpopsMeshRefinement/ to Path']);
addpath([currdir,'/lib/gpopsMeshRefinement/'],'-begin');

% add Common/
disp(['Adding Directory ',currdir,'/lib/gpopsCommon/ to Path']);
addpath([currdir,'/lib/gpopsCommon/'],'-begin');

% add gpopsUtilities/
disp(['Adding Directory ',currdir,'/gpopsUtilities/ to Path']);
addpath([currdir,'/gpopsUtilities/'],'-begin');

% add license/
disp(['Adding Directory ',currdir,'/license/ to Path']);
addpath([currdir,'/license/'],'-begin');

% add NLP solver directory
% $$$ if isdir('nlp/snopt'),
% $$$   disp(['Adding Directory ',currdir,'/nlp/snopt/ to Path']);
% $$$   addpath([currdir,'/nlp/snopt/'],'-begin');
% $$$ end
if isdir('nlp/snopt'),
  disp(['Adding Directory ',currdir,'/nlp/snopt/snopt2015/ to Path']);
  addpath([currdir,'/nlp/snopt/snopt2015'],'-begin');
end
if isdir('nlp/ipopt'),
  disp(['Adding Directory ',currdir,'/nlp/ipopt/ to Path']);
  addpath([currdir,'/nlp/ipopt/'],'-begin');
end

disp('-------------------------------------------------------------------------');
disp('The GPOPS-II directories have been successfully added to the MATLAB path.');
disp('-------------------------------------------------------------------------');
disp('');
pathNotSaved = savepath;
if pathNotSaved;
  warning('% The MATLAB path could not be saved to the master path definition file PATHDEF.M.    %');
  warning('% In order to include the GPOPS-II directories automatically each time MATLAB starts, %');
  warning('% please see the instructions in the GPOPS-II user guide.                             %');
end;

disp('');

if exist('adigator.m','file') < 2;
  disp('%--------------------------------------------------------------------------%');
  disp('% In order to install ADiGator, please visit the AdiGator project website: %')
  disp('% <a href="http://sourceforge.net/projects/adigator/">http://sourceforge.net/projects/adigator/</a>.                               %')
  disp('%--------------------------------------------------------------------------%');
end
  cd ./adigator
  startupadigator
  cd ..



