function [ind_op,ind_cp] = find_parenthesis_indices(str)

% Find all parenthesis indices (used to determine ending of cos/sin
% arguments, etc)
k_op = strfind(str,'('); % open parenthesis
k_cp = strfind(str,')'); % close parenthesis

% Determine indices of matching open and close parenthesis
k_p_sort = sort([k_op,k_cp]);

ind_op = nan(1,length(k_op));
ind_cp = nan(1,length(k_cp));
ctr_op = 0; % Array counter
for ctr = 1 : 1 : length(k_p_sort)
    
    if ismember(k_p_sort(ctr),k_op)
        % Add to op index
        ctr_op = ctr_op + 1;
        ind_op(ctr_op) = k_p_sort(ctr);
    else
        % Determine corresponding array location of closing parenthesis. Start at
        % end of open parenthesis array and work backwards until find last open
        % parenthesis.
        I = find(isnan(ind_op) == 0);
        while_ctr = 0;
        while 1
            if isnan(ind_cp(I(end-while_ctr)))
                ind_cp(I(end-while_ctr)) = k_p_sort(ctr);
                break;
            else
                while_ctr = while_ctr + 1;
            end
        end
    end
    
end

return
