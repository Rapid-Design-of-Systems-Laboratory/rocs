function [ out_str ] = makeSTTrhsKernels( eom_func, in, varargin)
% Generates C code for the EOM file from the symbolic expressions 
%
%   Takes state equation expressions and also generates expressions for 
%   n-th order Jacobians (depending on config setting)
%
% eom_func   : Valid function handle to system dynamics file
% config     : Structure containing configuration information
%
% out_str    : Contains the final C code output
%
    syms x t;
    out_str = '';
    
    % Copy in required variables from configuration
    auxdata_in  = in.const;
    nStates     = eom_func();
    maxSttOrder = in.gpuSolve.sttOrder;

    % Create symbolic variables for the states
    x = sym('x',[nStates 1]);
    
    % Create symbolic variables for the auxiliary data
    auxdata_fields = fieldnames(auxdata_in);
    auxdata_sym = struct();
    for i=1:numel(auxdata_fields)
        fieldname = char(auxdata_fields(i));
        auxdata_sym.(fieldname) = sym(fieldname);
    end
    in.const_sym = auxdata_sym;
    in.const_fields = auxdata_fields;
    
    % Evaluate the state equations
    fx = eom_func(t,x,auxdata_sym, varargin{:});
	
	out_str = sprintf('%sextern "C" __global__ void stateKernel(double *d_data, double *d_dx)\n{\n',out_str);
    % Write the auxiliary data
    for i=1:numel(auxdata_fields)
        fieldname = char(auxdata_fields(i));
        c_fieldname = sym2cexp(auxdata_sym.(fieldname),true);
        out_str = sprintf('%s\tconst double %s = %f;\n',out_str,c_fieldname,auxdata_in.(fieldname));
    end
    out_str = sprintf('%s\n',out_str);
    % Initialize the state variables
    % (Uses zero-based indexing)
    for i=1:nStates
        out_str = sprintf('%s\tdouble x%d = d_data[%d];\n',out_str,i,i-1);
    end
    out_str = strcat(out_str,['\n\tint idx = threadIdx.x + blockIdx.x * blockDim.x;\n']);
	
	% Write out the state equations
	out_str = sprintf('%s\tswitch(idx)\n\t{\n',out_str);
	for i=1:nStates
	    cexp = sym2cexp(fx(i));
		out_str = sprintf('%s\t\tcase %d:\n\t\t\td_dx[idx] = %s\n\t\t\tbreak;\n',out_str,i-1,cexp);
	end
	out_str = strcat(out_str,'\t}\n');
	out_str = strcat(out_str,'\n}\n');
	
    % Write out the partial derivative (Jacobian) expressions
    Fstar = fx;
    offset = nStates;
    for m=1:maxSttOrder
		fprintf('STT Order %d...\n',m);
        % Calculate nth order jacobians
        Fstar = jacobian(Fstar,x);
        out_str = sprintf('%sextern "C" __global__ void stt%d_F_kernel(double *d_data, double *d_F)\n{\n',out_str,m);
        % Write the auxiliary data
        for i=1:numel(auxdata_fields)
            fieldname = char(auxdata_fields(i));
            c_fieldname = sym2cexp(auxdata_sym.(fieldname),true);
            out_str = sprintf('%s\tconst double %s = %f;\n',out_str,c_fieldname,auxdata_in.(fieldname));
        end
        out_str = sprintf('%s\n',out_str);

        % Initialize the state variables
        % (Uses zero-based indexing)
        for i=1:nStates
            out_str = sprintf('%s\tdouble x%d = d_data[%d];\n',out_str,i,i-1);
        end

        if m == 1
            out_str = strcat(out_str,[...
                '\n\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n', ...
                '\tint j = threadIdx.y + blockIdx.y * blockDim.y;\n\n']);
				
            out_str = sprintf('%s\tif(i < %d && j < %d)\n',out_str,nStates,nStates);
            out_str = strcat(out_str,['\t{\n', ...
                    '\t\t// Special case for 1st order derivative\n', ...
                    '\t\t// Get 1-D index to get the derivative\n', ...
                    '\t\t// Column major form\n']);
            out_str = sprintf('%s\t\tint offset = %d + i + j*%d;\n\t\td_F += i+ j*%d;\n\n',out_str,nStates,nStates,nStates);
        else
            out_str = strcat(out_str,[...
                '\n\tint idx = threadIdx.x + blockIdx.x * blockDim.x;\n']);
			out_str = sprintf('%s\tif(idx < pow((double)%d,%d+1.0))\n',out_str,nStates,m);
			out_str = strcat(out_str,[...
					'\t{\n', ...
                    '\t\t// Get 1-D index to get the derivative\n']);
			out_str = sprintf('%s\t\tint offset = %d*(pow((double)%d,(double)%d)-1)/(%d-1) + idx;\n\t\td_F += idx;\n\n',out_str,nStates,nStates,m,nStates);
        end
        out_str = strcat(out_str,'\n\t\tswitch(offset)\n\t\t{\n');
        
        zero_terms = [];
        for i=1:numel(Fstar)
            dim_Fstar = nStates*ones(m+1,1);
            cexp = sym2cexp(Fstar(ind2sub(dim_Fstar,i)));
            if(strcmp(cexp,'0.0;'))
                zero_terms = [zero_terms offset+i-1];
            else
                out_str = sprintf('%s\t\t\tcase %d:\n\t\t\t\t*d_F = %s\n\t\t\t\tbreak;\n',out_str,offset+i-1,cexp);
            end
        end
        for ind = zero_terms
            out_str = sprintf('%s\t\t\tcase %d:\n',out_str,ind);
        end
        out_str = sprintf('%s\t\t\t\t*d_F = 0.0;\n\t\t\t\tbreak;\n',out_str);
        offset = offset + nStates^(m+1);
        out_str = sprintf('%s\t\t}\n\t}\n}\n',out_str);
    end
end

