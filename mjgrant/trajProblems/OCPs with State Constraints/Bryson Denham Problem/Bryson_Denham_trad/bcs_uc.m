% This is the boundary conditions file for Constrained Bryson Denham Problem
% input : Normalized time and state parameters of the trajectory
% output : Rate of change of state parameters with time
% Developed by : Kshitij Mall
% Last modified: 8 Feb, 2016

function res = bcs_uc(YL,YR, p)

% Calculate control at the end
u1 = 0;
ham1_end = YR(4,end)*YR(2,end) + YR(5,end)*u1^3 + 0.5*u1^2 + YR(6,end);

u2 = -1/(3*YR(5,end));
ham2_end = YR(4,end)*YR(2,end) + YR(5,end)*u2^3 + 0.5*u2^2 + YR(6,end);

if ham1_end<ham2_end
    ham_end = ham1_end;
else
    ham_end = ham2_end;
end

% u_end = -YR(5,end);

% Calculate hamitonian values at the end
% ham_end = YR(4,end)*YR(2,end) + YR(5,end)*u_end + 0.5*u_end^2 + YR(6,end);

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Residue %%
%%%%%%%%%%%%%%%%%%%%%%%
res = [ YL(1,1) - 0 % h(0) = 0
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
    ham_end] % Since the final time is free (not really as we made t a state as well), hamiltonian is 0

return

% End of file
