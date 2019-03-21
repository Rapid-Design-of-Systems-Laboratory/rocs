% This is the boundary conditions file for Constrained Bryson Denham Problem
% input : Normalized time and state parameters of the trajectory
% output : Rate of change of state parameters with time
% Developed by : Kshitij Mall
% Last modified: 18 Mar, 2016

function res = bcs_jumps(YL,YR, p)
% Calculate control at different points
u1_minus = 1;
u1_plus = -2*YL(5,2);
u2_minus = -2*YR(5,2);
u2_plus = -1;
u3_minus = -1;
u3_plus = -2*YL(5,end);
u_end = -2*YR(5,end);

% Calculate hamitonian values at different points of the multiarcs problem
ham1_minus = YR(4,1)*YR(2,1) + YR(5,1)*(-YR(1,1) + YR(2,1)*(1.4 - 0.14*YR(2,1)^2) + 4*u1_minus) + YR(6,1) + u1_minus^2 + YR(1,1)^2;
ham1_plus = YL(4,2)*YL(2,2) + YL(5,2)*(-YL(1,2) + YL(2,2)*(1.4 - 0.14*YL(2,2)^2) + 4*u1_plus) + YL(6,2) + u1_plus^2 + YL(1,2)^2;
ham2_minus = YR(4,2)*YR(2,2) + YR(5,2)*(-YR(1,2) + YR(2,2)*(1.4 - 0.14*YR(2,2)^2) + 4*u2_minus) + YR(6,2) + u2_minus^2 + YR(1,2)^2;
ham2_plus = YL(4,3)*YL(2,3) + YL(5,3)*(-YL(1,3) + YL(2,3)*(1.4 - 0.14*YL(2,3)^2) + 4*u2_plus) + YL(6,3) + u2_plus^2 + YL(1,3)^2;
ham3_minus = YR(4,3)*YR(2,3) + YR(5,3)*(-YR(1,3) + YR(2,3)*(1.4 - 0.14*YR(2,3)^2) + 4*u3_minus) + YR(6,3) + u3_minus^2 + YR(1,3)^2;
ham3_plus = YL(4,end)*YL(2,end) + YL(5,end)*(-YL(1,end) + YL(2,end)*(1.4 - 0.14*YL(2,end)^2) + 4*u3_plus) + YL(6,end) + u3_plus^2 + YL(1,end)^2;
ham_end = YR(4,end)*YR(2,end) + YR(5,end)*(-YR(1,end) + YR(2,end)*(1.4 - 0.14*YR(2,end)^2) + 4*u_end) + YR(6,end) + u_end^2 + YR(1,end)^2;

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Residue %%
%%%%%%%%%%%%%%%%%%%%%%%
res = [YL(1,1) + 5 % x1(0) = -5
    YL(2,1) + 5  % x2(0) = -5
    YL(3,1) - 0 % t(0) = 0
    YL(4,1) - p(5) % lamx1(0) = nu_x1_0
    YL(5,1) - p(6) % lamx2(0) = nu_x2_0
    YL(6,1) - p(7) % lamt(0) = nu_t_0
    YR(1,end) - p(8) % h(1) = nu_x1_f
    YR(2,end) - p(9) % v(1) = nu_x2_f
    YR(3,end) - 4.5 % t(1) = 4.5
    YR(4,end) - 0 % lamh(1) = nu_h_f
    YR(5,end) - 0 % lamv(1) = nu_v_f
    YR(6,end) - p(10) % lamt(1) = nu_t_f
    ham1_minus-ham1_plus % Hamiltonian does not jump between arc1 and arc2
    ham2_minus-ham2_plus % Hamiltonian does not jump between arc2 and arc3
    ham3_minus-ham3_plus % Hamiltonian does not jump between arc4 and arc4
    ham_end % Since the final time is free (not really as we made t a state as well), hamiltonian is 0
    YR(1,1) - YL(1,2) % Continuity of x1 at entry
    YR(2,1) - YL(2,2) % Continuity of x2 at entry
    YR(3,1) - YL(3,2) % Continuity of time at entry
    YR(4,1) - YL(4,2) % Continuity of lamx1 at entry
    %YR(5,1) + 0.5 % Continuity of lamx2 at entry
    YR(5,1) - YL(5,2) % Continuity of lamt at entry 
    YR(6,1) - YL(6,2) % Continuity of lamt at entry 
    YR(1,2) - YL(1,3) % Continuity of x1 at exit
    YR(2,2) - YL(2,3) % Continuity of x2 at exit  
    YR(3,2) - YL(3,3) % Continuity of time at exit 
    YR(4,2) - YL(4,3) % Continuity of lamx1 at exit
    YR(5,2) - YL(5,3) % Continuity of lamx2 at exit
    %YR(5,2) - 0.5
    YR(6,2) - YL(6,3) % Continuity of lamt at exit
    YR(1,3) - YL(1,4) % Continuity of x1 at exit
    YR(2,3) - YL(2,4) % Continuity of x2 at exit  
    YR(3,3) - YL(3,4) % Continuity of time at exit 
    YR(4,3) - YL(4,4) % Continuity of lamx1 at exit
    %YR(5,3) - 0.5
    YR(5,3) - YL(5,4) % Continuity of lamx2 at exit 
    YR(6,3) - YL(6,4)
    ]; % Continuity of lamt at exit 

return

% End of file
