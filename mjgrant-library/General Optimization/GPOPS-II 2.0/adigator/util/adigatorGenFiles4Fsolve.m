function funcs = adigatorGenFiles4Fsolve(setup)
% funcs = adigatorGenFiles4Fsolve(setup)
%
% ADiGator Jacobian File Generation Function: this
% function is used when you wish to generate derivative files for use in
% Matlab's optimization toolbox. The user must specify their objective
% function. This function calls adigator on the user's functions and
% creates a objective Hessian file of the form required by matlabs
% optimization toolbox. Note: all auxiliary data must be both known and
% fixed when passed to this function.
%
% --------------------------- User Function Requirements ---------------- %
% The objective function must be of the form:
%                   f = myfun(x) OR f = myfun(x,auxdata)
% Where x corresponds to the decision variable, and auxdata is a fixed
% structure/cell/array.
%
% If auxdata is an input, then the auxdata must be given to
% adigatorGenFiles4Fsolve.
%
% --------------------------- Input Structure --------------------------- %
%
% setup should be a structure with the following fields:
%
% setup.numvar:     length of the variable of differentiation (i.e. x)
%
% setup.function:  string name of the user's constraint function
%
% setup.options:    (optional) options structure generated by
%                   adigatorOptions
%
% setup.auxdata:    (optional) if the constraint function has an auxdata
%                   input, then the auxdata should be included here
%
%
% --------------------------- Output Structure -------------------------- %
%
% funcs.jacobian: jacobian function handle
%
% Jacobian function has output [fun, grd, hes]
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% see also adigator, adigatorCreateDerivInput, adigatorOptions,
% adigatorGenJacFile, adigatorGenHesFile, adigatorGenFiles4Fmincon,
% adigatorGenFiles4Fminunc

% ---------------------------- Parse Inputs ----------------------------- %

if isfield(setup,'numvar') && length(setup.numvar) == 1
  n = setup.numvar;
  if n < 2
    error('adigatorGenFiles4Fsolve not coded for scalar decision variable')
  end
else
  error('must specify decision vector length')
end
if isfield(setup,'function') && ischar(setup.function) && ...
    exist(setup.function,'file')
  ConFunName = setup.function;
  ConFun     = str2func(ConFunName);
else
  error('must specify function')
end
if isfield(setup,'options')
  opts = setup.options;
  if ~isfield(opts,'overwrite')
    opts.overwrite = 1;
  end
else
  opts.overwrite = 1;
end
if isfield(setup,'auxdata')
  if nargin(ConFun) ~= 2
    error('if auxdata specified, constraint function must have 2 inputs')
  end
  auxflag = 1;
  auxdata = setup.auxdata;
else
  if nargin(ConFun) ~= 1
    error('constraint function should have single input')
  end
  auxflag = 0;
end
if nargout(ConFun)~=1
  error('constraint function should have single output')
end

% --------------------------- Set Up File Names ------------------------- %
ConD1FileName  = [ConFunName,'_ADiGatorJac'];    % Con 1st derivs
JacFileName    = [ConFunName,'_Jac'];             % Jacobian Wrapper
AllFileNames   = {ConD1FileName, JacFileName};


% ---------------------- Check/Delete All Files ------------------------- %
CallingDir = cd;
for I = 1:length(AllFileNames)
  FileNamei = AllFileNames{I};
  if exist([CallingDir,filesep,FileNamei,'.m'],'file');
    if opts.overwrite
      delete([CallingDir,filesep,FileNamei,'.m']);
      rehash
    else
      error(['The file ',CallingDir,filesep,FileNamei,'.m already exists, ',...
        'quitting transformation. To set manual overwrite of file use ',...
        '''''adigatorOptions(''OVERWRITE'',1);''''. Alternatively, delete the ',...
        'existing file and any associated .mat file.']);
    end
  end
end

% ----------------------- Differentiate Constraint File ----------------- %
x = adigatorCreateDerivInput([n 1],'x');
if auxflag
  Inputs = {x, auxdata};
else
  Inputs = {x};
end
conout1 = adigator(ConFunName,Inputs,ConD1FileName,opts);
conout1i = conout1{1};
mi = prod(conout1i.func.size);

% -------------------------- Create Necessary Files --------------------- %
if auxflag
  InVarStr  = 'x,auxdata';
  dInVarStr = 'gx,auxdata';
else
  InVarStr  = 'x';
  dInVarStr = 'gx';
end

Jfid = fopen([JacFileName,'.m'],'w+');
fprintf(Jfid,['function [f, J] = ',JacFileName,'(',InVarStr,')\n']);

% Print Function Header
fprintf(Jfid,'%% \n');
fprintf(Jfid,'%% Wrapper file generated by ADiGator\n');
fprintf(Jfid,['%% ',char(169),'2010-2014 Matthew J. Weinstein and Anil V. Rao\n']);
fprintf(Jfid,'%% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ \n');
fprintf(Jfid,'%% Contact: mweinstein@ufl.edu\n');
fprintf(Jfid,'%% Bugs/suggestions may be reported to the sourceforge forums\n');
fprintf(Jfid,'%%                    DISCLAIMER\n');
fprintf(Jfid,'%% ADiGator is a general-purpose software distributed under the GNU General\n');
fprintf(Jfid,'%% Public License version 3.0. While the software is distributed with the\n');
fprintf(Jfid,'%% hope that it will be useful, both the software and generated code are\n');
fprintf(Jfid,'%% provided ''AS IS'' with NO WARRANTIES OF ANY KIND and no merchantability\n');
fprintf(Jfid,'%% or fitness for any purpose or application.\n\n');

% ---------------------------- Function Calls  -------------------------- %
cD1 = conout1i.deriv.nzlocs;
cD1nnz = size(cD1,1); 
fprintf(Jfid,'if nargout == 1\n');
fprintf(Jfid,['    f = ',ConFunName,'(',InVarStr,');\n']);
fprintf(Jfid,'elseif nargout == 2\n');
fprintf(Jfid,'    gx.f = x;\n');
fprintf(Jfid,'    gx.dx = ones(%1.0d,1);\n',n);
fprintf(Jfid,['   con = ',ConD1FileName,'(',dInVarStr,');\n']);
fprintf(Jfid,'    f = con.f;\n');
if mi == 0
  fprintf(Jfid,'    J = [];\n');
elseif cD1nnz == mi*n
  fprintf(Jfid,'    J = reshape(con.dx,%1.0f,%1.0f);\n',mi,n);
elseif mi == 1
  fprintf(Jfid,'    J = zeros(%1.0f,%1.0f);\n',mi,n);
  fprintf(Jfid,'    J(con.dx_location(:,1)) = con.dx;\n');
elseif mi*n >= 250 && cD1nnz/(mi*n) < 3/4
  fprintf(Jfid,['    J = sparse(con.dx_location(:,1),con.dx_location(:,2),',...
    'con.dx,%1.0f,%1.0f);\n'],mi,n);
else
  fprintf(Jfid,'    J = zeros(%1.0f,%1.0f);\n',mi,n);
  fprintf(Jfid,'    J((con.dx_location(:,2)-1)*%1.0f+con.dx_location(:,1)) = con.dx;\n',mi);
end
fprintf(Jfid,'end\n\n');

% --------------------------- Close All Files --------------------------- %
fclose(Jfid);
rehash

% -------------------------- Create Function Calls ---------------------- %
if auxflag
  funcs.jacobian = eval(['@(x)',JacFileName,'(x,auxdata)']);
else
  funcs.jacobian = str2func(JacFileName);
end

fprintf('\n<strong>adigatorGenFiles4Fsolve</strong> successfully generated Jac wrapper file\n\n');
