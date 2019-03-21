function [] = writeConstants(fid,const,writeIn,fileformat)

% Obtain fieldnames of constants
names = fieldnames(const);

if nargin < 4
	fileformat = 'matlab';
end

% Write constants to file
switch fileformat
	case 'matlab'
		fprintf(fid,'%% Constants\n');
		for ctr = 1 : 1 : length(names)
			if writeIn
				fprintf(fid,[char(names{ctr}),' = in.const.',char(names{ctr}),';\n']);
			else
				fprintf(fid,[char(names{ctr}),' = const.',char(names{ctr}),';\n']);
			end
		end
	case 'c'
		fprintf(fid,'// Constants\n');
		for ctr = 1 : 1 : length(names)
			fprintf(fid,['\tconst double ',char(names{ctr}),' = d_const[',num2str(ctr-1),'];\n']);
		end
end

return
