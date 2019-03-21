function [] = writeVehicleParameters(fid,veh,writeOut)

% Get input data
temp = veh.Cl;
varData = whos('temp');

% Obtain fieldnames of constants
names = fieldnames(veh);

% Write constants to file
fprintf(fid,'%% Vehicle Parameters\n');
for ctr = length(names) : -1 : 1 % Backwards so Aref written before Cl and Cd
    if strcmp(varData.class,'sym')
        fprintf(fid,[char(names{ctr}),' = ',char(veh.(names{ctr})),';\n']); % Make sure nonzero and always positive. Get transition from slender to blunt body. Need to add vehicle constraints later.
    else
        if writeOut
            fprintf(fid,[char(names{ctr}),' = out.VEH.',char(names{ctr}),';\n']);
        else
            fprintf(fid,[char(names{ctr}),' = VEH.',char(names{ctr}),';\n']);
        end
    end
end

return