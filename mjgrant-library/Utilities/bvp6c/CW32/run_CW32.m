function run_CW32(tol)
clc

% which tolerances to use
for tol = [1e-3]
    fid = fopen(['bvp6c.solve.errors' num2str(tol) '.txt'],'w');
    fprintf(fid,'prob,E,L2_Y,L2_DY,L_inf_Y,L_inf_DY,time,ODE evals,NumPoints,ignored\n');
    fprintf(fid,'tol,%2.2g\n',tol);
    fclose(fid);
    
    fid = fopen(['bvp5c.solve.errors' num2str(tol) '.txt'],'w');
    fprintf(fid,'prob,E,L2_Y,L2_DY,L_inf_Y,L_inf_DY,time,ODE evals,NumPoints,ignored\n');
    fprintf(fid,'tol,%2.2g\n',tol);
    fclose(fid);

    fid = fopen(['bvp4c.solve.errors' num2str(tol) '.txt'],'w');
    fprintf(fid,'prob,E,L2_Y,L2_DY,L_inf_Y,L_inf_DY,time,ODE evals,NumPoints,ignored\n');
    fprintf(fid,'tol,%2.2g\n',tol);
    fclose(fid);

    % which problems to solve
    for j=[1:32]
        CW32(j,tol);
        fprintf(' \n');
    end
end

