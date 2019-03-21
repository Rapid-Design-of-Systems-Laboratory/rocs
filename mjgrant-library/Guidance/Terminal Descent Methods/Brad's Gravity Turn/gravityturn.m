function [m_prop, m_final, T_required, exit_flag] = gravityturn(X0, Isp, CL, CD, S, alt_target,vel_target,Tmax,omega_mag,r_planet,mu_planet)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %DESCRIPTION:
% %
% %Determines the thrust to achieve a target altitude given an initial state 
% %and full 3DOF EOMs in an inertial frame to account for rotation and a 
% %spherical planet
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %INPUTS:
% % X0: Initial Inertial State Vector [x; y; z; xdot; ydot; zdot; mass]=[r; v] [m; m/s]
% % Isp: Specific Impulse [s]
% % CL: Coefficient of Lift [-]
% % CD: Coefficient of Drag [-]
% % S: Reference Area [m^2]
% % alt_target: Target Altitude [m]
% % vel_target: Relative Velocity Magnitude [m/s]
% % Tmax: Maximum Thrust Level [N]
% % omega_mag: Planetary Rotation Rate [rad/s]
% % r_planet: Planet Surface Radius [m]
% % mu_planet: Gravitational Parameter of the Planet [kg-m^3/s^2]
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %OUTPUTS:
% %
% % m_prop: Propellant mass [kg]
% % m_final: Final vehicle mass [kg]
% % T_required: Thrust required to achieve desired altitude at desired
% %             velocity [N]
% % exit_flag: 1 = Nominal; 0 = Not Converged
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%TEST CASE&&&&&&&&
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events',@events);
tspan = [0 5000];
%%%%%%%%%%%%%%%%%%%%%%%

%User Adjustable Parameters
alt_tolerance = 0.1;    %Convergence Altitude Tolerance
MaxIter = 1000;         %Maximum Number of Iterations for the Newton Iteration
timestep = 0.01;        %Timestep for Integration

%Compute Initial Thrust Guess
alt0 = norm(X0(1:3,1)) - r_planet;
vel0 = X0(4:6,1) - cross(omega_mag*[0;0;1],X0(1:3,1));
vel0_mag = norm(vel0);
FPA0 = acos(norm(cross(X0(1:3,1),X0(4:6,1)))/(norm(X0(1:3,1))*norm(X0(4:6,1))));

if dot(X0(1:3,1),X0(4:6,1)) < 0 && FPA0 > 0
    FPA0 = -FPA0;
elseif dot(X0(1:3,1),X0(4:6,1)) > 0 && FPA0 < 0
    FPA0 = -FPA0;
end

R0 = -alt0/sin(FPA0);
Tguess = X0(7,1)*(3.71*vel0_mag^2/(2*R0));

%Initialize Parameters
exit_flag = 1;
T = Tguess;
alt_final = alt_tolerance + 1000;
iter = 0;
omega = omega_mag*[0;0;1];

while abs(alt_final-alt_target) > alt_tolerance && iter < MaxIter
    
    %Determine final altitude for given thrust
    [t,X,te,Xe,ie] = ode45(@gravityturn_propagation,tspan,X0,options,T,Isp,CL,CD,S,vel_target,omega_mag,r_planet,mu_planet);

    te
    
    if ~isempty(Xe)
        Xe = Xe';
    else
        Xe = X(end,:)';
        te = NaN;
    end
    
    r_ii_final = Xe(1:3,1);
    r_ii_mag_final = norm(r_ii_final);
    alt_final = r_ii_mag_final - r_planet;
    
    v_ii_final = Xe(4:6,1);
    v_ii_mag_final = norm(v_ii_final);
    v_pp_final = v_ii_final - cross(omega,r_ii_final);
    v_pp_mag_final = norm(v_pp_final);
    
    m_final = Xe(7,1);
    m_prop = X0(7,1) - Xe(7,1);
    
    %Store variables for tracking
    alt_track(iter+1,1) = alt_final;
    T_track(iter+1,1) = T;
    m_prop_track(iter+1,1) = m_prop;
    time_track(iter+1,1) = te;
    
    if iter == 0
        
        T = 1.01*T; %For the first time through multiply the thrust by 1%
    
    else

        %Determine itermidiate vatriables
        T_new = T_track(iter+1,1);
        T_old = T_track(iter,1);
        alt_new = alt_track(iter+1,1);
        alt_old = alt_track(iter,1);
        
        if T_new == T_old
            exit_flag = 0;
            break;
        end
       
        %Check to see if a change in altitude occurred, if not decrease
        %step size, if so, perform a Newton step
        if alt_old == alt_new && timestep > 1e-16
            timestep = timestep/2;
        elseif timestep < 1e-16
            exit_flag = 0;
            break;
        else
            dT_dalt = (T_new-T_old)/(alt_new-alt_old);
            eta = alt_target - alt_new;
            T = T + dT_dalt*eta;
        end
    
        %If thrust is less than zero, reduce thrust by a factor of 2
        if T < 0
            T = Tguess/2;
            Tguess = Tguess/2;
        end
        
        if T > Tmax
            T = Tmax * 0.99;
            Tguess = Tmax * 0.99;
        end
             
    end

    iter = iter + 1;
end

if (iter == MaxIter || exit_flag == 0)
    exit_flag = 0;
    T_required = NaN;
    m_prop = NaN;
    m_final = NaN;
else
    exit_flag = 1;
    T_required = T;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dXdt] = gravityturn_propagation(t,X,T,Isp,CL,CD,S,vel_target,omega_mag,r_planet,mu_planet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs and Initial Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[force_ii] = traj_force(t,X,T,Isp,CL,CD,S,vel_target,omega_mag,r_planet,mu_planet);

% Derivatives
dXdt(1:3,1) = X(4:6,1);
dXdt(4:6,1) = force_ii/X(7,1);
dXdt(7,1) = -T/(9.81*Isp);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction] = events(t,X,T,Isp,CL,CD,S,vel_target,omega_mag,r_planet,mu_planet)

omega = omega_mag*[0;0;1];
r_ii = X(1:3,1);
v_ii = X(4:6,1);
v_ii_mag = norm(v_ii);
v_pp = v_ii - cross(omega,r_ii);
v_pp_mag = norm(v_pp);

value = v_pp_mag - vel_target;
isterminal = 1;
direction = -1;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [force_ii] = traj_force(t,X,T,Isp,CL,CD,S,vel_target,omega_mag,r_planet,mu_planet)
	
  pos_ii = X(1:3,1);                    % Inertial position
  pos_ii_mag = norm(pos_ii);            % Inertial position magnitude
  vel_ii = X(4:6,1);                    % Inertial velocity
  vel_ii_mag = norm(vel_ii);            % Inertial velocity magnitude
  mass = X(7,1);                        % Mass
  
  alt = pos_ii_mag - r_planet;          % Altitude
  [rho,pres] = get_rho_pres(alt);       % Density
 
  % Inertial to planet relative transformation
  rot_angle = omega_mag*t;                  % Rotation angle [rad]
  L_PI = [ cos(rot_angle) sin(rot_angle) 0; ...
          -sin(rot_angle) cos(rot_angle) 0; ...
                        0              0 1];
  
  pos_ii_hat = pos_ii/pos_ii_mag;       % Inertial position vector direction
  h_ii = cross(pos_ii,vel_ii);          % Inertial angular momenum vector [m^2/s]
  h_ii_mag = norm(h_ii);                % Magnitude of inertial angular momenum vector [m^2/s]
  pos_pp = L_PI*pos_ii;                 % Position vector planet/planet [m]
  omega = omega_mag * [0;0;1];
  vel_pp = L_PI*(vel_ii - cross(omega,pos_ii)); % Velocity vector planet/planet [m/s]
  vel_pp_mag = norm(vel_pp);
  vel_pp_hat = vel_pp/norm(vel_pp);
  h_pp = cross(pos_pp,vel_pp);
  h_pp_mag = norm(h_pp);
  h_pp_hat = h_pp/h_pp_mag;
    
  % Dynamic pressure
  q = 1/2*rho*vel_pp_mag^2;

  % Obtain guidance values
  bank = 0;                                 % Bank angle [deg]
  thrust_pp = -T*(vel_pp/vel_pp_mag);
  thrust_ii = L_PI'*thrust_pp;              % Inertial thrust vector [N]

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Vectors and Rotation Tensors of Interset %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  n1 = skew(-h_pp_hat);
  R1 = eye(3) + sind(90)*n1 + (1-cosd(90))*n1*n1;
  n2 = skew(vel_pp_hat);
  R2 = eye(3) + sind(bank)*n2 + (1-cosd(bank))*n2*n2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Vehicle Aerodynamic Force Calculations %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  drag_pp_hat = -vel_pp_hat;                % Planet relative drag force direction
  drag_pp = q*CD*S*drag_pp_hat;             % Planet relative drag force vector
  lift_pp_hat = R2*R1*vel_pp_hat;           % Planet relative lift force direction
  lift_pp = q*CL*S*lift_pp_hat;             % Planet relative lift force vector
  drag_ii = L_PI'*drag_pp;                  % Inertial drag force vector
  lift_ii = L_PI'*lift_pp;                  % Inertial lift force vector
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Gravity Force Calculation %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gravity_ii_mag = mu_planet*mass/pos_ii_mag^2;
  gravity_ii = gravity_ii_mag*(-pos_ii_hat); % Inertial gravity force [N]

  %%%%%%%%%%%%%%%%%
  %% Total Force %%
  %%%%%%%%%%%%%%%%%

  % Total inertial external force vector on body [N]
  force_ii = drag_ii + lift_ii + gravity_ii + thrust_ii;
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho,press] = get_rho_pres(alt)

[rho_table,press_table] = get_atm();

rho = interp1q(rho_table(:,1),rho_table(:,2),alt);
press = interp1q(press_table(:,1),press_table(:,2),alt);

if isnan(rho)
    rho = rho_table(1,2);
end

if isnan(press)
    press = press_table(1,2);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n = skew(vec)

n = [     0  -vec(3)  vec(2); ...
      vec(3)      0  -vec(1); ...
     -vec(2)  vec(1)      0];

return
