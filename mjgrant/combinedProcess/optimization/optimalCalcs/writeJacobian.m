function writeJacobian(fid,varString,func,jacY)

fprintf(fid,[varString,' = [']);
for ctr = 1 : 1 : length(func)
    if ctr ~= 1
        % Add spaces
        fprintf(fid,'                   ');
    end
    % Write Jacobian
    jacSet = jacobian(func(ctr),jacY);
    numCol = length(jacSet);
    for col = 1 : 1 : numCol
        fprintf(fid,'%s',char(jacSet(col)));
        if col ~= numCol
            fprintf(fid,', ');
        end
    end
    % Close out line
    if ctr == length(func)
        % Last entry
        fprintf(fid,'];\n');
    else
        % Middle entry
        fprintf(fid,'; ...\n');
    end
end

return