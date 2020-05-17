function [] = writeQuantities(fid,quantity,fileformat)

% Obtain fieldnames of quantities
names = fieldnames(quantity);

if nargin < 3
	fileformat = 'matlab';
end

% Write quantities to file
switch fileformat
	case 'matlab'
		fprintf(fid,'%% Derived Quantities\n');
		for ctr = 1 : 1 : length(names)
			fprintf(fid,[char(names{ctr}),' = ',char(quantity.(names{ctr}){1}),';\n']);
		end
	case 'c'
		fprintf(fid,'// Derived Quantities\n');
		for ctr = 1 : 1 : length(names)
			fprintf(fid,['\tconst double ',char(names{ctr}),' = d_quantity[',num2str(ctr-1),'];\n']);
		end
end

return
