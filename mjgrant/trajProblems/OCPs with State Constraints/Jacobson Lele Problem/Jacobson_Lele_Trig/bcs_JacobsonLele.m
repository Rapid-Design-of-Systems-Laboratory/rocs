% This is the boundary conditions file for Constrained Bryson Denham Problem
% input : Normalized time and state parameters of the trajectory
% output : Rate of change of state parameters with time
% Developed by : Kshitij Mall
% Last modified: 8 Mar, 2016

function res = bcs_JacobsonLele(YL,YR, p)

% Calculate control at different points
u_init = -YR(9,1);
u_end1 = -YL(9,end);
u_end2 = -YR(9,end);

% Calculate hamitonian values at different points of the multiarcs problem
ham1_minus = YR(6,1)*YR(2,1) + YR(7,1)*YR(3,1) + YR(8,1)*YR(4,1) + YR(9,1)*u_init + 0.5*u_init^2 + YR(10,1);
ham2_plus = YL(6,2)*YL(2,2) + YL(7,2)*YL(3,2) + YL(8,2)*YL(4,2) + YL(9,2)*0 + 0.5*0^2 + YL(10,2);
ham3_minus = YR(6,2)*YR(2,2) + YR(7,2)*YR(3,2) + YR(8,2)*YR(4,2) + YR(9,2)*0 + 0.5*0^2 + YR(10,2);
ham4_plus = YL(6,end)*YL(2,end) + YL(7,end)*YL(3,end) + YL(8,end)*YL(4,end) + YL(9,end)*u_end1 + 0.5*u_end1^2 + YL(10,end);
ham_end = YR(6,end)*YR(2,end) + YR(7,end)*YR(3,end) + YR(8,end)*YR(4,end) + YR(9,end)*u_end2 + 0.5*u_end2^2 + YR(10,end);

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Residue %%
%%%%%%%%%%%%%%%%%%%%%%%
res = [ YL(1,1) - 0 % x1(0) = 0
    YL(2,1) - 15/12  % x2(0) = 15/12
    YL(3,1) + 15/12  % x3(0) = -15/12
    YL(4,1) - 15/16  % x4(0) = 15/16
    YL(5,1) - 0 % t(0) = 0
    YL(6,1) - p(4) % lamx1(0) = nu_x1_0
    YL(7,1) - p(5) % lamx2(0) = nu_x2_0
    YL(8,1) - p(6) % lamx3(0) = nu_x3_0
    YL(9,1) - p(7) % lamx4(0) = nu_x4_0
    YL(10,1) - p(8) % lamt(0) = nu_t_0
    YR(1,end) - 0 % x1(10) = 0
    YR(2,end) + 15/12 % x2(10) = -15/12
    YR(3,end) + 15/12 % x3(10) = -15/12
    YR(4,end) + 15/16 % x4(10) = 15/16
    YR(5,end) - 10 % t(10) = 10
    YR(6,end) - p(9) % lamx1(1) = nu_x1_f
    YR(7,end) - p(10) % lamx2(10) = nu_x2_f
    YR(8,end) - p(11) % lamx3(10) = nu_x3_f
    YR(9,end) - p(12) % lamx4(10) = nu_x4_f
    YR(10,end) - p(13) % lamt(10) = nu_t_f
    ham1_minus-ham2_plus % Hamiltonian does not jump between arc1 and arc2
    ham3_minus-ham4_plus % Hamiltonian does not jump between arc2 and arc3
    ham_end % Since the final time is free (not really as we made t a state as well), hamiltonian is 0
    YR(1,1) - YL(1,2) % Continuity of h at entry
    YR(2,1) - YL(2,2) % Continuity of v at entry
    YR(3,1) - YL(3,2) % Continuity of time at entry
    YR(4,1) - YL(4,2) % Continuity of time at entry
    YR(5,1) - YL(5,2) % Continuity of time at entry
    YR(1,2) - YL(1,3) % Continuity of h at exit
    YR(2,2) - YL(2,3) % Continuity of v at exit  
    YR(3,2) - YL(3,3) % Continuity of time at exit 
    YR(4,2) - YL(4,3) % Continuity of time at exit 
    YR(5,2) - YL(5,3) % Continuity of time at exit 
    YR(6,1) - YL(6,2) - p(14)% Discontinuity of lamh at entry
    YR(7,1) - YL(7,2) - p(15)% Discontinuity of lamv at entry
    YR(8,1) - YL(8,2) - p(16)% Dicontinuity of lamt at entry Note: No jump
    YR(9,1) - YL(9,2) - p(17)% Dicontinuity of lamt at entry Note: No jump
    YR(10,1) - YL(10,2) 
    YR(6,2) - YL(6,3) - 0% Continuity of lamh at exit
    YR(7,2) - YL(7,3) - 0% Continuity of lamv at exit 
    YR(8,2) - YL(8,3) - 0% Continuity of lamt at exit 
    YR(9,2) - YL(9,3) - 0% Continuity of lamt at exit 
    YR(10,2) - YL(10,3) - 0% Continuity of lamt at exit 
    YL(1,2) - 1 % Constraint Derivative 1 tangency condition
    YL(2,2) - 0
    YL(3,2) - 0
    YL(4,2) - 0];

return

% End of file
