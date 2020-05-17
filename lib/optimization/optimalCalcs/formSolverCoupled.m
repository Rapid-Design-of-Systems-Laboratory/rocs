function [solution, isNumericalSoln, numSolFun] = formSolverCoupled(expression,value,variable,assumptions,in,oc)
    % This function creates an expression that represents the solution to 'expression = value' for variable.
    % The supplied assumptions are used to try and construct an analytic solution. If no analytic solution
    % can be found, a numerical root solving (fzero) expression is constructed.
	
    % initialGuess = 0; % make this smarter in the future
    isNumericalSoln = false;
    numSolFun = '';
    
    for ctrControl = 1 : 1 : oc.num.controls
        % Determine if an analytic solution can be found
        sol = mathematica_solve(expression(ctrControl),sym('0'),variable(ctrControl),assumptions);

        % Check first solution (in case multiple solutions exist for Mathematica inverse function of unknown, user-supplied function)
        if isempty(strfind(char(sol(1)),'InverseFunction')) && isempty(strfind(char(sol(1)),'Solve')) % Analytic solution found, finish		
            % if length(char(solution(1))) > 1
            % solution not singular
            solution{ctrControl} = sol;
            % else % solution singular
            % solution{ctrControl} = solution;
            % end		
        else % Analytic solution not found, form numerical root solving expression
            isNumericalSoln = true;
            solution{ctrControl} = [['u(',num2str(ctrControl),')']];
            break;  % Proceed to numerical solution
        end
    end

    if isNumericalSoln
        % solution{1} = ['fzero(@(',char(variable),') ',char(expression),'-',char(value),',',num2str(initialGuess),')'];
		
        % Determine range for interval search. Start by searching for a function that uses variable as input.
        indFunc = 0;
        indInput = 0;
		
        stopSearch = false;
		
        if isfield(in,'customFunc')

        else
      
            interval = [-48*pi/180 48*pi/180]; %%% NEED TO GENERALIZE
      
        end
		
        % Write first part of solver line
        numSolFun = ['u = newtonRaphson(@controlSolutionFunction,'];
        numSolFun = [numSolFun,'ones(',num2str(oc.num.controls),',1)*',num2str(mean(interval)),',1e-6'];
        
        % solution{1} = ['fzero(@',char(variable),'SolutionFunction,0*pi/180,',in.fzeroOptions];
        % solution{1} = ['zeroByLooping(@',char(variable),'SolutionFunction,[',num2str(interval),']'];
		
        % Write remaining custom function to be passed into the function
        extraInput = [];
        % for ctrInput = 2 : 1 : length(customFunc(indFunc).input)
        % 	extraInput = [extraInput,',',char(customFunc(indFunc).input(ctrInput){1})];
        % end
		
        % Write general variables to pass to root solving function
        % 		extraInput = [extraInput,',xAndLambda,const,constraint,bank']; %%%%%%% FIX HARDCODED CONTROL!!!!!!!!!!!!!
        extraInput = [extraInput,',xAndLambda,const,constraint'];
		
        % Write full solution expression
        numSolFun = [numSolFun,extraInput,');'];
		
        % Write function that fzero will call
        fid = fopen([in.autocodeDirTraj,'/controlSolutionFunction.m'],'w');
		
        % Header
        fprintf(fid,['function [zeroVal] = controlSolutionFunction(u',extraInput,')\n']);
        % fprintf(fid,'coder.extrinsic(''keyboard'');\n'); %%%%%%%%%%%%%% DELETE LATER!!!!!!!!!!!!!!!!!!!!!!
        % Assign input values
        if isfield(in,'const')
            writeConstants(fid,in.const,false);
            fprintf(fid,'\n');
        end
	
        if isfield(in,'constraintVal')
            writeConstraints(fid,in.constraintVal,false);
        end
        fprintf(fid,'\n');
	
        % Use real values since these will be fed into the user custom functions. Some Matlab algorithms require real numbers
        % (e.g., interpolation routines).
        fprintf(fid,'\n');
        fprintf(fid,'%% States\n');
        for ctrState = 1 : 1 : oc.num.states
            fprintf(fid,[char(oc.state.var(ctrState)),' = real(xAndLambda(',int2str(ctrState),',1));\n']);
        end
        fprintf(fid,'\n');

        fprintf(fid,'%% Costates\n');
        for ctrCostate = 1 : 1 : oc.num.costates
            fprintf(fid,[char(oc.costate.var(ctrCostate)),' = real(xAndLambda(',int2str(ctrCostate+oc.num.states),',1));\n']);
        end
        fprintf(fid,'\n');
	
        fprintf(fid,'%% Controls\n');
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,[char(oc.control.var(ctrControl)),' = real(u(',int2str(ctrControl),',1));\n']);
        end
        fprintf(fid,'\n');
        
        % Write derived quantities
        if isfield(in,'quantity')
            writeQuantities(fid,in.quantity);
            fprintf(fid,'\n');
        end
		
        % Write root-solving expression
        fprintf(fid,['zeroVal = real([\n']);
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,[char(expression(ctrControl)),'-',char(value),'; ...\n']);
        end
        fprintf(fid,']);\n'); % real expression required by fzero
		
        % fprintf(fid,'if zeroVal == inf\n');%abs(real(zeroVal)) > 1e8 || abs(imag(zeroVal)) > 0\n');
        % fprintf(fid,'  xAndLambda\n');
        % fprintf(fid,'  alfa\n');
        % fprintf(fid,'  bank\n');
        % fprintf(fid,'  zeroVal = zeroVal\n');
        % 
        % % fprintf(fid,'  alfaSet = linspace(-10*pi/180,10*pi/180,100);\n');
        % % fprintf(fid,'  zeroValSet = NaN(size(alfaSet));\n');
        % % fprintf(fid,'  for ctr = 1 : 1 : length(alfaSet)\n');
        % % fprintf(fid,'    alfa = alfaSet(ctr);\n');
        % % fprintf(fid,['    zeroValSet(ctr) = ',char(expression),'-',char(value),';\n']);
        % % fprintf(fid,'  end\n');
        % % fprintf(fid,'plot(alfaSet,zeroValSet);\n');
        % % 
        % % fprintf(fid,'  keyboard\n');
        % fprintf(fid,'end\n');
		
        % fprintf(fid,'if abs(real(zeroVal)) < 1e-3\n');
        % fprintf(fid,'  alfa\n');
        % fprintf(fid,'  bank\n');
        % fprintf(fid,'  zeroVal = zeroVal\n');
        % fprintf(fid,'  keyboard\n');
        % fprintf(fid,'end\n');
		
        fprintf(fid,'\n');
        fprintf(fid,'return\n');
        fprintf(fid,'\n');

        % Close file
        fclose(fid);
		
    end
    % end
	
    return

