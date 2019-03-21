% This is the boundary conditions file for Constrained Bryson Denham Problem
% input : Normalized time and state parameters of the trajectory
% output : Rate of change of state parameters with time
% Developed by : Kshitij Mall
% Last modified: 18 Mar, 2016

function res = bcs_new(YL,YR, p)
% Calculate control at different points
u_end1 = -1;
ham1 = YR(4,end)*YR(2,end) + YR(5,end)*(-YR(1,end) + YR(2,end)*(1.4 - 0.14*YR(2,end)^2) + 4*u_end1) + YR(6,end) + u_end1^2 + YR(1,end)^2;

u_end2 = 1;
ham2 = YR(4,end)*YR(2,end) + YR(5,end)*(-YR(1,end) + YR(2,end)*(1.4 - 0.14*YR(2,end)^2) + 4*u_end2) + YR(6,end) + u_end2^2 + YR(1,end)^2;

if ham1<ham2
    ham_end = ham1;
else
    ham_end = ham2;
end
u_end3 = -2*YR(5,end);
ham3 = YR(4,end)*YR(2,end) + YR(5,end)*(-YR(1,end) + YR(2,end)*(1.4 - 0.14*YR(2,end)^2) + 4*u_end3) + YR(6,end) + u_end3^2 + YR(1,end)^2;

if ham3<ham_end
    ham_end = ham3;
end

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Residue %%
%%%%%%%%%%%%%%%%%%%%%%%
res = [YL(1,1) + 5 % x1(0) = -5
    YL(2,1) + 5  % x2(0) = -5
    YL(3,1) - 0 % t(0) = 0
    YL(4,1) - p(2) % lamx1(0) = nu_x1_0
    YL(5,1) - p(3) % lamx2(0) = nu_x2_0
    YL(6,1) - p(4) % lamt(0) = nu_t_0
    YR(1,end) - p(5) % h(1) = nu_x1_f
    YR(2,end) - p(6) % v(1) = nu_x2_f
    YR(3,end) - 4.5 % t(1) = 4.5
    YR(4,end) - 0 % lamh(1) = nu_h_f
    YR(5,end) - 0 % lamv(1) = nu_v_f
    YR(6,end) - p(7) % lamt(1) = nu_t_f
    ham_end]; % Since the final time is free (not really as we made t a state as well), hamiltonian is 0

return

% End of file