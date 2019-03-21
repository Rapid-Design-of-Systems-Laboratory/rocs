function [value,isterminal,direction] = gravity_turn_event_func(t,y,param)

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

h = y(3) - 1; % Altitude - end at 10 m above ground

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign Event Values %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Terminate at h = 0
value = h; % set altitude as variable of interest
isterminal = 1; % set integration to terminate
direction = -1; % h should be descending

return

