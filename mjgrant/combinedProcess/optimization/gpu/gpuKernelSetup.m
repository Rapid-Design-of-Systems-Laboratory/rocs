function [ptxfile] = gpuKernelSetup(eom_func, in, out, varargin)  
	% Generate the C code for the given EOM function
	% Returns the absolute path to the CUDA ptx assembly file generated
   
	nStatesAndCostates = 2*in.oc.num.states;
	
	fname = strcat(pwd,'/',func2str(eom_func),'_stt',num2str(in.gpuSolve.sttOrder));
	fext = in.gpuSolve.eomFormat;
	
	eomfile = strcat(fname,'.',fext);
	% if in.gpuSolve.generateEoms ~= true
	% 	if exist(strcat(fname,'.',fext),'file')
	% 		ptxfile = strcat(fname,'.',fext);
	% 		return;
	% 	end
	% end
	fprintf('\nGenerating equations for GPU...\n');
	if in.gpuSolve.sttOrder == 1
		out_str = makeSplitStmKernel(eom_func,in,out,varargin{:});
	else
		out_str = makeSTTrhsKernels(eom_func,in,varargin{:});
	end
	fprintf('Creating CU file...\n');	
	cufile = strcat(fname,'.cu');
	
	% Open the output file
	fp = fopen(cufile,'w+');
	fprintf(fp,out_str);
	fclose(fp);
	
	% Generate assert statements for all constants
	const_str = '';
	auxdataFields = fieldnames(in.const);
	for i=1:length(auxdataFields)
		const_str = sprintf('%s\tassert(isa(constants.%s,''double''));\n',const_str,char(auxdataFields(i)));
	end
	
	if isfield(in,'constraintVal')
		constraint_names = fieldnames(in.constraintVal);
	else
		constraint_names = {};
	end

	% Write constants to file
	const_str = sprintf('%s\tassert(isa(constraint,''struct''));\n',const_str);
	const_str = strcat(const_str,'');

	for ctr = 1 : 1 : length(constraint_names)
		const_str = sprintf('%s\tassert(isa(constraint.%s,''double''));\n',const_str,char(constraint_names{ctr}));
	end
	
	% Output split EOM propagation function
	spliteomfile = strcat(pwd,'/',func2str(eom_func),'_split.m');
	out_str = fileread('splitstate_func.tmpl.m');
	out_str = strrep(out_str,'__NSTATES__',num2str(nStatesAndCostates));
	out_str = strrep(out_str,'__MAXARCS__',num2str(in.maxNumArcs));
	
	out_str = strrep(out_str,'__CONSTANTS__',const_str);
	out_str = strrep(out_str,'__EOMNAME__',func2str(eom_func));
	fp = fopen(spliteomfile,'w');
	fwrite(fp,out_str);
	fclose(fp);
	% End of spilt EOM function
	
	% Output RK4 EOM propagation function
	rk4eomfile = strcat(pwd,'/',func2str(eom_func),'_RK4.m');
	out_str = fileread('stateRK4_func.tmpl.m');
	out_str = strrep(out_str,'__NSTATES__',num2str(nStatesAndCostates));
	out_str = strrep(out_str,'__MAXARCS__',num2str(in.maxNumArcs));
	out_str = strrep(out_str,'__CONSTANTS__',const_str);
	out_str = strrep(out_str,'__EOMNAME__',func2str(eom_func));
	fp = fopen(rk4eomfile,'w');
	fwrite(fp,out_str);
	fclose(fp);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Run codegen on the MATLAB RK4 files
	codegen(spliteomfile);
	codegen(rk4eomfile)
	
	% Compile the NVIDIA CUDA files
	curpath = pwd;

	nvcccmd = sprintf('nvcc %%s -m64 -ccbin g++ %s -o %%s -Xptxas -O0 -Xcompiler -Wno-unused-variable -I %s',cufile,fileparts(mfilename('fullpath')));
	if strcmp(in.gpuSolve.eomFormat,'cubin')
		ptxfile = strcat(fname,'.cubin');
		nvcccmd = sprintf(nvcccmd,'-cubin -arch=sm_20',ptxfile);
	else
		ptxfile = strcat(fname,'.ptx');
		nvcccmd = sprintf(nvcccmd,'-ptx',ptxfile);
	end
	nvcccmd
	if in.gpuSolve.write
		system(nvcccmd);
	end
	%delete(cufile);
	
 %	if ~exist(strcat(fileparts(mfilename('fullpath')),'propagator.',mexext()),'file') || ...
%			~exist(strcat(fileparts(mfilename('fullpath')),'propagator_stm_gpu.',mexext()),'file') 
%		copyfile(strcat(fileparts(mfilename('fullpath')),'/../build/propagator.',mexext()),'.');
%		copyfile(strcat(fileparts(mfilename('fullpath')),'/../build/propagator_stm_gpu.',mexext()),'.')
%	end
end