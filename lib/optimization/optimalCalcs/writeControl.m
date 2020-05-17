function [] = writeControl(in,oc,control,coefficients,funcName,hamiltonian,controlWriteOrder)

% Open file to write control function. Also used in postprocessing to determine optimal control history.
fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
% fid_c = fopen([in.autocodeDirTraj,'/',funcName,'_.h'],'w'); % This line causes an error with Matlab's codegen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Header

% if strcmpi(in.gpuSolve.derivativeMethod,'csd')
% fprintf(fid_c,'#include <cusp/complex.h>\n#define complex_t cusp::complex<double>\n#include<stdio.h>\n');
% end

% fprintf(fid_c,['__device__ _num_t ',funcName,'Hamiltonian(_num_t *xAndLambda,double *d_const,double *d_constraints,']);

fprintf(fid,'function [');
for ctrControl = 1 : 1 : oc.num.controls
	fprintf(fid,[char(oc.control.var(ctrControl)),'Save']);
	% fprintf(fid_c,['_num_t ',char(oc.control.var(ctrControl))]);
	if ctrControl ~= oc.num.controls
		fprintf(fid,',');
		% fprintf(fid_c,',');
	end
end
% fprintf(fid_c,[')\n{\n']);
fprintf(fid,[',hamiltonianSave] = ',funcName,'(xAndLambda,const,constraint,numArcs)\n']);

fprintf(fid,'coder.extrinsic(''keyboard'');\n');

% fprintf(fid,['if ~isa(xAndLambda,''sym'')\n']);
writeHeader(fid,in);
fprintf(fid,'assert(isa(xAndLambda,''double''));\n');
if ~in.convertParametersToStates
	fprintf(fid,['assert(all(size(xAndLambda)== [',int2str(oc.num.states*2),',1]));\n']);
else
	fprintf(fid,['assert(all(size(xAndLambda)== [',int2str(oc.num.states*2+in.maxNumArcs),',1]));\n']);
	% fprintf(fid,['coder.varsize(''xAndLambda'',[',int2str(oc.num.states*2+in.maxNumArcs),',1],[1,0]);\n']);
end
% fprintf(fid,['end\n\n']);
fprintf(fid,'assert(isa(numArcs,''double''));\n');
fprintf(fid,['assert(all(size(numArcs)== [1,1]));\n']);

% Assign input values
if isfield(in,'const')
	writeConstants(fid,in.const,false);
	% writeConstants(fid_c,in.const,false,'c');
end
fprintf(fid,'\n\n');
% fprintf(fid_c,'\n\n');
if isfield(in,'constraintVal')
	writeConstraints(fid,in.constraintVal,false);
	% writeConstraints(fid_c,in.constraintVal,false,'c');
end
fprintf(fid,'\n\n');
% fprintf(fid_c,'\n\n');

fprintf(fid,'%% States\n');
% fprintf(fid_c,'\t// States\n');
for ctrState = 1 : 1 : oc.num.states
	fprintf(fid,[char(oc.state.var(ctrState)),' = complex(xAndLambda(',int2str(ctrState),',1));\n']);
	% fprintf(fid_c,['\t_num_t ',char(oc.state.var(ctrState)),' = xAndLambda[',int2str(ctrState-1),'];\n']);
end
fprintf(fid,'\n');

fprintf(fid,'%% Costates\n');
for ctrCostate = 1 : 1 : oc.num.costates
	fprintf(fid,[char(oc.costate.var(ctrCostate)),' = complex(xAndLambda(',int2str(ctrCostate+oc.num.states),',1));\n']);
	% fprintf(fid_c,['\t_num_t ',char(oc.costate.var(ctrCostate)),' = xAndLambda[',int2str(ctrCostate+oc.num.states-1),'];\n']);
end
fprintf(fid,'\n');
% fprintf(fid_c,'\n');

% if in.convertParametersToStates
% 	fprintf(fid,'%% Independant Variables\n');
% 	fprintf(fid,['tSet = NaN(',int2str(in.maxNumArcs),',1);\n']);
% 	fprintf(fid,['for ctr = 1 : 1 : numArcs\n']);
% 	fprintf(fid,['\ttSet(ctr) = xAndLambda(',int2str(2*oc.num.states),'+ctr);\n']);
% 	fprintf(fid,['end\n']);
% 	fprintf(fid,'\n');
% 	
% 	fprintf(fid_c,['\t_num_t t_[',int2str(in.maxNumArcs),'];\n']);
% 	fprintf(fid_c,'\tfor(int ctr=0;ctr<blockDim.z;ctr++)\n\t{\n');
% 	fprintf(fid_c,['\t\tt_[ctr] = xAndLambda[',int2str(2*oc.num.states),'+ctr];\n']);
% 	fprintf(fid_c,'\t}\n');
% end

% if in.rootSolving == 0
	fprintf(fid,'%% Control\n');
	for ctrControl = 1 : 1 : oc.num.controls
		fprintf(fid,[char(oc.control.var(ctrControl)),' = NaN;\n']);
		fprintf(fid,[char(oc.control.var(ctrControl)),'Save = NaN;\n']);
		% fprintf(fid_c,['_num_t ',char(oc.control.var(ctrControl)),' = ',char(oc.control.var(ctrControl)),'_.real();']);
	end
	fprintf(fid,'\n');
	
	% Write derived quantities
	if isfield(in,'quantity')
		writeQuantities(fid,in.quantity);
		fprintf(fid,'\n');
	end
% end

% fprintf(fid_c,['\t_num_t hamiltonian = ',sym2cexp(hamiltonian,true),';\n']);
% fprintf(fid_c,['\t_num_t hamiltonian = ',sym2cexp(hamiltonian,true),';\n']);
% fprintf(fid_c,['\t_num_t hamiltonian = 0;\n']);
% fprintf(fid_c,['printf("%%lf\\n",hamiltonian);\n']);
% fprintf(fid_c,'\treturn hamiltonian;\n');


% fprintf(fid_c,'\n}\n');

% Compute control value and determine sufficiency condition
% fprintf(fid,['ctrlSave = NaN(',int2str(oc.numControl),',1) + 0*i;\n']);
% [m,n] = size(Huu);
% 	fprintf(fid,['d2Hdu2 = NaN(',int2str(m),',',int2str(n),');\n']);


% fprintf(fid_c,'template <typename T>\n');
% fprintf(fid_c,['__device__ void ',funcName,'(_num_t *xAndLambda,double *d_const,double *d_constraints,']);

for ctrControl = 1 : 1 : oc.num.controls
	% fprintf(fid_c,['_num_t *',char(oc.control.var(ctrControl)),'Save,']);
end
% fprintf(fid_c,'_num_t *hamiltonianSave)\n{\n');

% writeConstants(fid_c,in.const,false,'c');
% fprintf(fid_c,'\n\n');
if isfield(in,'constraintVal')
	% writeConstraints(fid_c,in.constraintVal,false,'c');
end
% fprintf(fid_c,'\n\n');

% fprintf(fid_c,'\t// States\n');
% for ctrState = 1 : 1 : oc.num.states
% 	fprintf(fid_c,['\t_num_t ',char(oc.state.var(ctrState)),' = xAndLambda[',int2str(ctrState-1),'];\n']);
% end
% fprintf(fid_c,'\t// Costates\n');
% for ctrCostate = 1 : 1 : oc.num.costates
% 	fprintf(fid_c,['\t_num_t ',char(oc.costate.var(ctrCostate)),' = xAndLambda[',int2str(ctrCostate+oc.num.states-1),'];\n']);
% end
% fprintf(fid_c,'\n');

% fprintf(fid_c,'\t// Control\n');
% for ctrControl = 1 : 1 : oc.num.controls
% 	fprintf(fid_c,['\t_num_t ',char(oc.control.var(ctrControl)),';\n']);
% 	% fprintf(fid_c,['\tdouble ',char(oc.control.var(ctrControl)),'Save;\n']);
% end
% fprintf(fid_c,'\n');

% Define/initialize Hamiltonian
if in.minimize
	fprintf(fid,'hamiltonian = inf;\n');
	fprintf(fid,'hamiltonianSave = inf;\n');
	% fprintf(fid_c,'\tdouble hamiltonianSave = INFINITY;\n double hamiltonian;\n');
	% fprintf(fid_c,'\t_num_t hamiltonian;\n');
else
	fprintf(fid,'hamiltonian = -inf;\n');
	fprintf(fid,'hamiltonianSave = -inf;\n');
	% fprintf(fid_c,'\tdouble hamiltonianSave = -INFINITY;\n double hamiltonian;\n');
	% fprintf(fid_c,'\t_num_t hamiltonian;\n');
end


% Determine number of possible controls and test each possible combination. Select best control based on Pontryagin's Minimum Principle.
numControlOptions = NaN(1,oc.num.controls);
for ctrControl = 1 : 1 : oc.num.controls
	numControlOptions(ctrControl) = length(control{ctrControl});
end
numControlCombinations = prod(numControlOptions);
controlIndex = cell(1,oc.num.controls);

fprintf(fid,'I = 1i;\n');

% % Determine maximum number of coefficients to declare in C source file
% maxCoeffCount = 0;
% for ctrControl = 1 : 1 : oc.num.controls
%     controlWriteIndex = controlWriteOrder(ctrControl);
%     for ctr = 1 : 1 : length(coefficients{controlWriteIndex})
%         for ctr2 = 1 : 1 : length(coefficients{controlWriteIndex}{ctr})
%             maxCoeffCount = max(maxCoeffCount, length(coefficients{controlWriteIndex}{ctr2}));
%         end
%     end
% end
% 
% if maxCoeffCount > 0
% 	fprintf(fid_c,['\tcomplex_t ']);
% 	for ctrCoefficient = 1 : 1 : maxCoeffCount
% 		fprintf(fid_c,['coeff',int2str(ctrCoefficient)]);
% 		if ctrCoefficient == maxCoeffCount
% 			fprintf(fid_c,';\n');
% 		else
% 			fprintf(fid_c,', ');
% 		end
% 	end
% end
% fprintf(fid_c,'\tcomplex_t I_(0,1);\n');
% 
% if oc.num.controls > 0
% 	fprintf(fid_c,['\tcomplex_t ',]);
% 	for ctrControl = 1 : 1 : oc.num.controls
% 		fprintf(fid_c,[char(oc.control.var(ctrControl)),'_']);
% 		if ctrControl == oc.num.controls
% 			fprintf(fid_c,';\n');
% 		else
% 			fprintf(fid_c,', ');
% 		end
% 	end
% end

% Write each control combination
for ctrControlCombinations = 1 : 1 : numControlCombinations
	% ctrControlCombinations = 4;
	
	% Determine indices of each control array to write
	[controlIndex{:}] = ind2sub(numControlOptions,ctrControlCombinations);
	
	for ctrControl = 1 : 1 : oc.num.controls
		
		controlWriteIndex = controlWriteOrder(ctrControl);
		% Compute coefficient values
		fprintf(fid,'%% Coefficients\n');
% 		fprintf(fid_c,'\t// Coefficients\n');
    for ctrCoefficient = 1 : 1 : length(coefficients{controlWriteIndex}{controlIndex{controlWriteIndex}})
			fprintf(fid,['coeff',int2str(ctrCoefficient),' = (', ...
 				char(coefficients{controlWriteIndex}{controlIndex{controlWriteIndex}}(ctrCoefficient)),');\n']);
% 			fprintf(fid_c,['\tcoeff',int2str(ctrCoefficient),' = ', ...
% 				sym2cexp(coefficients{controlWriteIndex}{controlIndex{controlWriteIndex}}(ctrCoefficient),true),';\n']);
    end
    
    
    
    
    
    
    
    
    if in.rootSolving == 0
      
      % Perform a check to see if control variables coupled. If they are, cannot
      % solve them in succession.
      coupledExpression = false;
%       for ctrControlCheck = 1 : 1 : oc.num.controls
%         
%         try
%           char(control{controlWriteIndex}(controlIndex{controlWriteIndex}))
%         catch
%           keyboard
%         end
%         
%         if strfind(char(oc.control.var(ctrControlCheck)), ...
%           char(control{controlWriteIndex}(controlIndex{controlWriteIndex})))
%           coupledExpression = true;
%         end
%       end

%% Change made in March 2015: If no control found then assign it a zero value     
      if ~coupledExpression
      if isempty(char(control{controlWriteIndex}(controlIndex{controlWriteIndex})))
          % fprintf(fid,[char(oc.control.var(controlWriteIndex)),' = solve(coeff1,',char(oc.control.var(controlWriteIndex)),',''Real'',true);\n']);
          fprintf(fid,[char(oc.control.var(controlWriteIndex)),' = 0;\n']);
        % disp('Oh no! No control!')
      else
        % Write expression
        fprintf(fid,[char(oc.control.var(controlWriteIndex)),' = real(',char(control{controlWriteIndex}(controlIndex{controlWriteIndex})),');\n']);
      end
      else
        % Numerically solve for roots
      end
    else
			fprintf(fid,[char(oc.control.var(controlWriteIndex)),' = (',char(control{controlWriteIndex}(controlIndex{ctrControl})),');\n']);
    end
    
%     try
%     fprintf(fid_c,['\t',char(oc.control.var(controlWriteIndex)),'_ = (',sym2cexp(sym(control{controlWriteIndex}(controlIndex{ctrControl})),true),');\n']);
% 		fprintf(fid_c,['\t',char(oc.control.var(controlWriteIndex)),' = ',char(oc.control.var(controlWriteIndex)),'_.real();\n']);
%     catch
%       fprintf('Need to fix c expression construction.\n');
%     end
    
    
%         fprintf(fid,'if abs(lamV) < 1e-8\n');
%         fprintf(fid,'  ang = 0;\n');
%         fprintf(fid,'else\n');
% %         fprintf(fid,'  ang = abs(ang);\n');
%         fprintf(fid,'end\n');
    
    
    
    
    
    
		% Write protection statement for singularities
		fprintf(fid,['if isnan(real(',char(oc.control.var(controlWriteIndex)),')) || isinf(real(',char(oc.control.var(controlWriteIndex)),'))\n']);
		% fprintf(fid_c,['\tif(isNaN(',char(oc.control.var(controlWriteIndex)),'))\n\t{\n']);
		fprintf(fid,['\t',char(oc.control.var(controlWriteIndex)),' = 0;\n']);
		% fprintf(fid_c,['\t\t',char(oc.control.var(controlWriteIndex)),' = 0;\n']);
		fprintf(fid,'end\n');
		% fprintf(fid_c,'\t}\n');
    


    
    
%     fprintf(fid,'if abs(alfa) > 100*pi/180\n');
%     fprintf(fid,'  alfa = 100*pi/180*sign(alfa);\n');
%     fprintf(fid,'end\n');
    
    
    
	end
	
	% Check Pontryagin's Minimum Principle to select proper control
		
	% Compute Hamiltonian
	% fprintf(fid_c,[hamiltonianStr,'\n']);
	fprintf(fid,['hamiltonian = real(',char(hamiltonian),');\n']);
	

% 	% fprintf(fid_c,['\thamiltonian = ',funcName,'Hamiltonian(xAndLambda,d_const,d_constraints,']);
% 	for ctrControl = 1 : 1 : oc.num.controls
% 		fprintf(fid_c,[char(oc.control.var(ctrControl))]);
% 		if ctrControl ~= oc.num.controls
% 			fprintf(fid_c,',');
% 		end
% 	end
% 	fprintf(fid_c,[');\n']);

	% fprintf(fid_c,['\thamiltonian = 0;\n']);
	% fprintf(fid_c,['\thamiltonian = ',sym2cexp(hamiltonian,true),';\n']);
	% Save controls if minimizing/maximizing Hamiltonian
	if ctrControlCombinations > 1
		if in.minimize
			fprintf(fid,'if hamiltonian < hamiltonianSave\n');
			% fprintf(fid_c,'\tif(compare(hamiltonian,*hamiltonianSave,1e-4)==-1)\n\t{\n');
		else
			fprintf(fid,'if hamiltonian > hamiltonianSave\n');
			% fprintf(fid_c,'\tif(compare(hamiltonian,*hamiltonianSave,1e-4)==1)\n\t{\n');
		end
	end
	for ctrControl = 1 : 1 : oc.num.controls
		fprintf(fid,['\t',char(oc.control.var(ctrControl)),'Save = real(',char(oc.control.var(ctrControl)),');\n']);
		% fprintf(fid_c,['\t\t*',char(oc.control.var(ctrControl)),'Save = ',char(oc.control.var(ctrControl)),';\n']);
	end
	% fprintf(fid,['disp(''Control #',num2str(ctrControlCombinations),''')\n']);
	% fprintf(fid_c,['printf("Control #',num2str(ctrControlCombinations),'\\n");\n']);
	
	fprintf(fid,'\thamiltonianSave = hamiltonian;\n');
	% fprintf(fid_c,'\t\t*hamiltonianSave = hamiltonian;\n');
	
	if ctrControlCombinations > 1
		fprintf(fid,'end\n');
		% fprintf(fid_c,'\t}\n');
	end
		
	% [m,n] = size(Huu);
	% 	
	% 			for ctr1 = 1 : 1 : m
	% 				for ctr2 = 1 : 1 : n
	% 		
	% 				% Write sufficiency condition. Do after compute control values since highly coupled. ...
	% 				% Compute for bank only for now
	% 				fprintf(fid,['d2Hdu2(',int2str(ctr1),',',int2str(ctr2),') = real(',char(Huu(ctr1,ctr2)),');\n']);
	% 	
	% 			end
	% end
	% else
	% fprintf(fid,['d2Hdu2(',int2str(ctr1),',1) = 0;\n']);	
end

% fprintf(fid,'ctrlSave1\n');
% fprintf(fid,'ctrlSave2\n');
% fprintf(fid,'d2Hdu2\n');

% fprintf(fid,'bank = min(abs(ctrlSave1));\n');
% fprintf(fid,[char(oc.u(2)),'=real(',char(oc.uOptimal.unconstrained{2}(1)),');\n']);

% fprintf(fid,'bankVal = min(abs(ctrlSave1))\n');
% fprintf(fid,'[Y,I] = max(d2Hdu2)\n');
% fprintf(fid,'I = 2;\n')
% fprintf(fid,'bank = ctrlSave1(I(1))\n');
% fprintf(fid,'bankSave = bank;\n')
% fprintf(fid,'alfa = ctrlSave2(I(1))\n');
% fprintf(fid,[char(oc.u(2)),'=real(',char(oc.uOptimal.unconstrained{2}(1)),')\n']);

% fprintf(fid,'if isnan(bank)\n');
% fprintf(fid,'  bank = 0;\n');
% fprintf(fid,['  ',char(oc.u(2)),'=real(',char(oc.uOptimal.unconstrained{2}(1)),');\n']);
% fprintf(fid,'end\n');
% fprintf(fid,'if isnan(alfa)\n');
% fprintf(fid,'  alfa = 0;\n');
% fprintf(fid,'end\n');

% fprintf(fid,'if bankVal ~= bankSave && ~isnan(bankVal)\n');
% fprintf(fid,'bankVal\n');
% fprintf(fid,'bankSave\n');
% fprintf(fid,'error(''stop'');\n');
% fprintf(fid,'end\n');

% fprintf(fid,'bank = 0;\n');
% fprintf(fid,[char(oc.u(2)),'=real(',char(oc.uOptimal.unconstrained{2}(1)),');\n']);

% fprintf(fid,'error(''stop here'')\n');

% fprintf(fid,['d2Hdu2 = NaN(',int2str(oc.numControl),',1);\n']);
% for ctr1 = 1 : 1 : length(oc.u)
% 	
% 	% Write sufficiency condition. Do after compute control values since highly coupled
% 	fprintf(fid,['d2Hdu2(',int2str(ctr1),',1) = real(',char(oc.Huu.unconstrained(ctr1)),');\n']);
% 	
% end

for ctrControl = 1 : 1 : oc.num.controls
	% fprintf(fid_c,['\t\tprintf("',char(oc.control.var(ctrControl)),' : %%lf\\n",',char(oc.control.var(ctrControl)),'Save.real());\n']);
	% fprintf(fid,[char(oc.control.var(ctrControl)),'Save\n']);
end

% fprintf(fid,'keyboard');

fprintf(fid,'\n');
fprintf(fid,'return\n');
fprintf(fid,'\n');

% fprintf(fid_c,'\n');
% fprintf(fid_c,'\tswitch(ctrl_index)\n\t{\n');
% for ctrControl = 1 : 1 : oc.num.controls
% 	fprintf(fid_c,['\t\tcase ',num2str(ctrControl-1),':\n']);
% 	fprintf(fid_c,['\t\t\treturn ',char(oc.control.var(ctrControl)),'Save;\n']);
% end
% fprintf(fid_c,'\t}\n');
% fprintf(fid_c,'\treturn 0;\n');
% fprintf(fid_c,'}\n');
% Close file
fclose(fid);
% fclose(fid_c);

return

