% This is the boundary conditions file for Constrained Bryson Denham Problem
% input : Normalized time and state parameters of the trajectory
% output : Rate of change of state parameters with time
% Developed by : Kshitij Mall
% Last modified: 8 Feb, 2016

function res = bcs_jumps(YL,YR, p)
% Calculate control at different points
u_init = -YR(5,1);
u_end1 = -YL(5,end);
u_end2 = -YR(5,end);

% Calculate hamitonian values at different points of the multiarcs problem
ham1_minus = YR(4,1)*YR(2,1) + YR(5,1)*u_init + 0.5*u_init^2 + YR(6,1);
ham2_plus = YL(4,2)*YL(2,2) + YL(5,2)*0 + 0.5*0^2 + YL(6,2);
ham3_minus = YR(4,2)*YR(2,2) + YR(5,2)*0 + 0.5*0^2 + YR(6,2);
ham4_plus = YL(4,3)*YL(2,3) + YL(5,3)*u_end1 + 0.5*u_end1^2 + YL(6,3);
ham_end = YR(4,3)*YR(2,3) + YR(5,3)*u_end2 + 0.5*u_end2^2 + YR(6,3);

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Residue %%
%%%%%%%%%%%%%%%%%%%%%%%
res = [ YL(1,1) - 0 % h(0) = 0
    YL(2,1) - 1  % v(0) = 1
    YL(3,1) - 0 % t(0) = 0
    YL(4,1) - p(4) % lamh(0) = nu_h_0
    YL(5,1) - p(5) % lamv(0) = nu_v_0
    YL(6,1) - p(6) % lamt(0) = nu_t_0
    YR(1,end) - 0 % h(1) = 0
    YR(2,end) + 1 % v(1) = -1
    YR(3,end) - 1 % t(1) = 1
    YR(4,end) - p(7) % lamh(1) = nu_h_f
    YR(5,end) - p(8) % lamv(1) = nu_v_f
    YR(6,end) - p(9) % lamt(1) = nu_t_f
    ham1_minus-ham2_plus % Hamiltonian does not jump between arc1 and arc2
    ham3_minus-ham4_plus % Hamiltonian does not jump between arc2 and arc3
    ham_end % Since the final time is free (not really as we made t a state as well), hamiltonian is 0
    YR(1,1) - YL(1,2) % Continuity of h at entry
    YR(2,1) - YL(2,2) % Continuity of v at entry
    YR(3,1) - YL(3,2) % Continuity of time at entry
    YR(1,2) - YL(1,3) % Continuity of h at exit
    YR(2,2) - YL(2,3) % Continuity of v at exit  
    YR(3,2) - YL(3,3) % Continuity of time at exit 
    YR(4,1) - YL(4,2) - p(10)% Discontinuity of lamh at entry
    % YL(5,2) - 0
    YR(5,1) - YL(5,2) - p(11)% Discontinuity of lamv at entry
    YR(6,1) - YL(6,2) - p(12)% Dicontinuity of lamt at entry Note: No jump
    YR(4,2) - YL(4,3) - 0% Continuity of lamh at exit
    % YR(5,2) - 0
    YR(5,2) - YL(5,3) - 0% Continuity of lamv at exit 
    YR(6,2) - YL(6,3) - 0% Continuity of lamt at exit 
    YL(1,2) - 1/8 % Constraint Derivative 1 tangency condition
    YL(2,2) - 0
    YL(3,2) - 0];
%    YL(3,2) - 0.375]
%     YR(1,2) - 1/8 % Constraint Derivative 1 tangency condition
%     YR(2,2) - 0
%     YR(3,2) - 0.625]
%     YR(1,2) - 1/8
%     YR(2,2) - 0]% % Constraint Derivative 2 tangency condition

return

% End of file
