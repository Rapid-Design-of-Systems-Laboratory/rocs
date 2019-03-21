function [ cexp ] = sym2cexp( symexp, trunc_semicolons, data_type )
    if nargin < 2
        trunc_semicolons = false;
    end
	if numel(symexp) == 1
    	cexp = strtrim(ccode(symexp));
	else
		fname = ['expr',int2str(rand*100),'.tmp'];
		ccode(symexp,'file',fname);
		fp = fopen(fname,'r');
		cexp = fscanf(fp,'%s');
		fclose(fp);
		% delete(fname);
		% delete('expr.tmp');
	end
	if nargin < 3
		data_type = 'double';
	end
	if nargin < 2
		trunc_semicolons = false;
	end
	if numel(symexp) == 1
	    if trunc_semicolons
	        cexp = cexp(6:end-1);
	    else
	        cexp = cexp(6:end);
	    end
	else
		repl = sprintf('\t%s $1\n',data_type);
		cexp = regexprep(cexp,'(t\d+\=[^;]+?;)',repl);
	end

	% Use constants of appropriate precision
	cexp = regexprep(cexp,'(pow\([^,]*,)',strcat(['$1(',data_type,')']));
	cexp = strrep(cexp,'sqrt(-1.0)','I_');
	% Substitute original expression for arccos() 
	% cexp = strrep(cexp,'3.141592653589793*-1.0+acos(','-acos(-');
	% cexp = strrep(cexp,'3.141592653589793-acos(','acos(-');
end

