% This is the boundary conditions file for Constrained Bryson Denham Problem
% input : Normalized time and state parameters of the trajectory
% output : Rate of change of state parameters with time
% Developed by : Kshitij Mall
% Last modified: 6 Nov, 2016

function res = bcs_trig(YL,YR, p)
global x1Max
% Calculate control at different points
u_end = -YR(5,end);

% Calculate hamitonian values at different points of the multiarcs problem
ham_end = YR(4,end)*YR(2,end)/(x1Max*cos(YR(1,end))) + YR(5,end)*u_end + 0.5*u_end^2 + YR(6,end);

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Residue %%
%%%%%%%%%%%%%%%%%%%%%%%
res = [ YL(1,1) - 0% 0 % h(0) = 0
    YL(2,1) - 1  % v(0) = 1
    YL(3,1) - 0 % t(0) = 0
    YL(4,1) - p(2) % lamh(0) = nu_h_0
    YL(5,1) - p(3) % lamv(0) = nu_v_0
    YL(6,1) - p(4) % lamt(0) = nu_t_0
    YR(1,end) - 0 % h(1) = 0
    YR(2,end) + 1 % v(1) = -1
    YR(3,end) - 1 % t(1) = 1
    YR(4,end) - p(5) % lamh(1) = nu_h_f
    YR(5,end) - p(6) % lamv(1) = nu_v_f
    YR(6,end) - p(7) % lamt(1) = nu_t_f
    ham_end]; % Since the final time is free (not really as we made t a state as well), hamiltonian is 0

return

% End of file
