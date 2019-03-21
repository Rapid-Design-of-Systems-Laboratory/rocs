function [] = writeScalingParameters(fid,scale,writeIn)

% Obtain fieldnames of constants
names = fieldnames(scale);

% Write constants to file
fprintf(fid,'%% Scaling parameters\n');
for ctr = 1 : 1 : length(names)
    if writeIn
        fprintf(fid,[char(names{ctr}),' = in.scale.',char(names{ctr}),';\n']);
    else
        fprintf(fid,[char(names{ctr}),' = scale.',char(names{ctr}),';\n']);
    end
end

return