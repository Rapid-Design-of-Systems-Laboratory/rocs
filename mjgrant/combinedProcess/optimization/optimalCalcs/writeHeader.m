function [] = writeHeader(fid,in)

% Write autocoding header
fprintf(fid,'%%#codegen\n');

% Write constants
% Obtain fieldnames of constants
if isfield(in,'const')
	names = fieldnames(in.const);
else
	names = {};
end

% Write constants to file
fprintf(fid,'assert(isa(const,''struct''));\n');
for ctr = 1 : 1 : length(names)
    fprintf(fid,['assert(isa(const.',char(names{ctr}),',''double''));\n']);
end

% Write constraints
% Obtain fieldnames of constraints
if isfield(in,'constraintVal')
	names = fieldnames(in.constraintVal);
else
	names = {};
end

% Write constants to file
fprintf(fid,'assert(isa(constraint,''struct''));\n');
for ctr = 1 : 1 : length(names)
    fprintf(fid,['assert(isa(constraint.',char(names{ctr}),',''double''));\n']);
end

% Write scaling parameters
% Obtain fieldnames of scaling parameters
%names = fieldnames(in.scale);

% Write constants to file
%fprintf(fid,'assert(isa(scale,''struct''));\n');
%for ctr = 1 : 1 : length(names)
%    fprintf(fid,['assert(isa(scale.',char(names{ctr}),',''double''));\n']);
%end

% Write vehicle parameters
% Obtain fieldnames of vehicle parameters
%names = fieldnames(veh);

% Write constants to file
%fprintf(fid,'assert(isa(VEH,''struct''));\n');
%for ctr = 1 : 1 : length(names)
%    fprintf(fid,['assert(isa(VEH.',char(names{ctr}),',''double''));\n']);
    % fprintf(fid,['assert(isa(VEH.',char(names{ctr}),',''sym''));\n']);
	%end

return
