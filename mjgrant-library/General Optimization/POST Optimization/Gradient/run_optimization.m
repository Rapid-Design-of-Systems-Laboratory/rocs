function run_optimization

% RUN_OPTIMIZATION optimizes the closed-loop performance with dispersions.
%   The objective of the optimization is to deploy at the highest altitude
%   while having a maximum range error <= a desired input value.

%%%%%%%%%%%%
%% INPUTS %%
%%%%%%%%%%%%

% Files and folder
param.reference_file = 'ref_fp2-212.mat';  % Reference trajectory file 
param.input_deck = 'nom';  % (Input deck name).inp
param.inp_file = 'nom.inp';  % Closed-loop file
% The following out.folder will store all of your output
% Also, be sure nothing is important in the directory! (It is wiped clean)
% The out.folder should be used as a temporary storage of your run files.
% Move them to your PC when the simulations are finished!
param.out_folder = '/home/vab/mjgrant/MSL/05-17a/output/temp/';

% Optimization parameters
param.range_err_tol = 20;  % Maximum range error tolerance [km]
num_bank_points = 6;   % Number of evenly-spaced points for bank profile
upper_vel_bank = 5000; % Upper bound for bank velocity profile
lower_vel_bank = 0;    % Lower bound for bank velocity profile
ub = 180*ones(1,num_bank_points); % Upper bound for bank profile points
lb = -180*ones(1,num_bank_points); % Lower bound for bank profile points
x0 = 70*ones(1,num_bank_points); % Initial guess for bank profile
options = []; % Optimization options (May want to change tolerances here)

% Stress cases
%           range(steep   shallow) alt(1 2)
stress.gam_set = [-14.2   -13.6        1 1]; % Initial gamma
stress.tau_set = [0.1     0.9          1 1]; % Initial dusttau
stress.lod_set = [0.8     1.2          1 1] ; % Initial L/D
stress.wgt_set = 2606.0*9.8066*[0.8 1.2 1 1];  % Initial weight [N] on Earth's surface (ballistic coefficient)
stress.lat_set = [-43.8917071]; % Initial latitude
stress.lon_set = [269.800747]; % Initial longitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET VARIABLES/PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine velocity points of bank profile
inc = (upper_vel_bank - lower_vel_bank)/num_bank_points;
param.vel_points = [lower_vel_bank : inc : upper_vel_bank];

%%%%%%%%%%%%%%%%%%%%%%
%% RUN OPTIMIZATION %%
%%%%%%%%%%%%%%%%%%%%%%

[x,fval,EXITFLAG,output,LAMBDA,grad,hessian]=...
    fmincon('entry_performance',x0,[],[],[],[],lb,ub,'entry_performance',options,param,stress);

%%%%%%%%%%%%
%% OUTPUT %%
%%%%%%%%%%%%

% Matlab's optimization performance
fprintf('\n\n');
fprintf('Optimization Results:\n');
fprintf('  Minimum parachute deploy altitude: %g\n',fval);
fprintf('  Bank Profile:\n');
for counter = 1 : 1 : length(x);
    fprintf('    Velocity [m/s] = %g, Bank Angle [deg] = %g\n',param.vel_points(counter),x(counter));
end
fprintf('  Altitude gradient: %g\n',grad);
fprintf('  Hessian:           %g\n',hessian);
fprintf('\n');
fprintf('Optimization Performance:\n');
fprintf('  Number of iterations:    %g\n',output.iterations);
fprintf('  Number of POST runs:     %g\n',output.funcCount);
fprintf('  Algorithm used:          %g\n',output.algorithm);
fprintf('  Number of CG iterations: %g\n',output.cgiterations);
fprintf('  First order optimality:  %g\n',output.firstorderopt);
fprintf('  Optimization message:    %s\n',output.message);

fprintf('\n');




