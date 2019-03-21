function [out1,out2] = entry_performance(x,param,stress)
%
% ENTRY_PERFORMANCE is used to run POST and quanitify the entry performance
%   based on range error and altitude at parachute deployment
%
% [out1, out2] = entry_performance(x,param,stress)
%
% Input: x -- bank angles at given velocity points
%        param -- structure containing all parameters necessary to run function
%        stress -- structure containing all parameters associated with the stress cases
%
% Output: out1 -- may be objective function value or constraint value based on
%                   what fmincon is currently doing
%         out2 -- may be empty or constratint value based on what fmincon is
%                   currently doing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE REFERENCE PROFILE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run POST open-loop using bank profile -- reference profile

% Use reference profile to generate gains
generate_gains([param.input_deck,'.inp']);

%%%%%%%%%%%%%%%%%%%%%%
%% RUN STRESS CASES %%
%%%%%%%%%%%%%%%%%%%%%%

% Assume the number of stress cases is equal to the length of the gamma set
for counter = 1 : 1 : length(stress.gam_set)
    
    % Set up stress case parameters
    param.lat_set = stress.lat_set(counter); % Initial latitude
    param.lon_set = stress.lon_set(counter); % Initial longitude
    param.gam_set = stress.gam_set(counter); % Initial gamma
    param.tau_set = stress.tau_set(counter); % Initial dusttau
    param.lod_set = stress.lod_set(counter); % Initial L/D
    param.wgt_set = stress.wgt_set(counter); % Initial weight [N] on Earth's surface (ballistic coefficient)
    
    % Run stress case
    run_closed_trajectory(param);
    
    % Extract alitude and range error at parachute deploy
    data_post = load([param.out_folder,'post/',param.input_deck,'1.inp']);
    data_apollo = load([param.out_folder,'apollo1.out']);
    alt_set(counter) = data.hgtareoid(end);  % Altitude at chute deploy
    range_err_set(counter) = data_apollo(end,11);  % Range left to fly at chute deploy
    
end

% Extract largest range error and lowest altitude at parachute deploy
alt = min(alt_set);
range_err = max(range_err_set);

%%%%%%%%%%%%
%% OUTPUT %%
%%%%%%%%%%%%

if nargout == 1  % fmincon checking minimization
    out1 = -alt;  % fmincon requires a minimization problem
elseif nargout == 2  % fmincon checking constraints
    out1 = [range_err - range_err_tol];  % range error must be within tolerance
    out2 = [];  % no equality constraint
end




