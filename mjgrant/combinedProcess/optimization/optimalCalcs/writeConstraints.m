function [] = writeConstraints(fid,constraint,writeIn,fileformat)

% Obtain fielnames of constants
names = fieldnames(constraint);

if nargin < 4
	fileformat = 'matlab';
end

switch(fileformat)
case 'matlab'
	% Write constants to file
	fprintf(fid,'%% Constraints\n');
	for ctr = 1 : 1 : length(names)
		if writeIn
			fprintf(fid,[char(names{ctr}),' = in.constraint.',char(names{ctr}),';\n']);
		else
			fprintf(fid,[char(names{ctr}),' = constraint.',char(names{ctr}),';\n']);
		end
	end
case 'c'
	fprintf(fid,'// Constraints\n');
	for ctr = 1 : 1 : length(names)
		fprintf(fid,['\tdouble ',char(names{ctr}),' = d_constraints[',num2str(ctr-1),'];\n']);
	end	
end
return
