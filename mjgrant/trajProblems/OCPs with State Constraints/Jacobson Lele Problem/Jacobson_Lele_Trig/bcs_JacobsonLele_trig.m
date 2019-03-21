% This is the boundary conditions file for Unconstrained Jacobson Lele Problem
% input : Normalized time and state parameters of the trajectory
% output : Rate of change of state parameters with time
% Developed by : Kshitij Mall
% Last modified: 6 Nov, 2016

function res = bcs_JacobsonLele_trig(YL,YR, p)
global scale x1Max
% Calculate control at different points
u_end = -YR(9,end);
% ham_end = YR(6,end)*YR(2,end) + YR(7,end)*YR(3,end) + YR(8,end)*YR(4,end) + YR(9,end)*u_end + 0.5*u_end^2 + YR(10,end);

ham_end = YR(6,end)*YR(2,end)*pi*(1+YR(1,end)^2)/(2*scale*x1Max) + YR(7,end)*YR(3,end) + YR(8,end)*YR(4,end) + YR(9,end)*u_end + 0.5*u_end^2 + YR(10,end);

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Residue %%
%%%%%%%%%%%%%%%%%%%%%%%
res = [ YL(1,1) - 0 % x1(0) = 0
    YL(2,1) - 15/12  % x2(0) = 15/12
    YL(3,1) + 15/12  % x3(0) = -15/12
    YL(4,1) - 15/16  % x4(0) = 15/16
    YL(5,1) - 0 % t(0) = 0
    YL(6,1) - p(2) % lamx1(0) = nu_x1_0
    YL(7,1) - p(3) % lamx2(0) = nu_x2_0
    YL(8,1) - p(4) % lamx3(0) = nu_x3_0
    YL(9,1) - p(5) % lamx4(0) = nu_x4_0
    YL(10,1) - p(6) % lamt(0) = nu_t_0
    YR(1,end) - 0 % x1(10) = 0
    YR(2,end) + 15/12 % x2(10) = -15/12
    YR(3,end) + 15/12 % x3(10) = -15/12
    YR(4,end) + 15/16 % x4(10) = 15/16
    YR(5,end) - 10 % t(10) = 10
    YR(6,end) - p(7) % lamx1(1) = nu_x1_f
    YR(7,end) - p(8) % lamx2(10) = nu_x2_f
    YR(8,end) - p(9) % lamx3(10) = nu_x3_f
    YR(9,end) - p(10) % lamx4(10) = nu_x4_f
    YR(10,end) - p(11) % lamt(10) = nu_t_f
    ham_end - 0]; % Since the final time is free (not really as we made t a state as well), hamiltonian is 0

return

% End of file
