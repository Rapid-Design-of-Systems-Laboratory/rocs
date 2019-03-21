function phaseout = rlvEntryContinuous(input)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
rad  = input.phase.state(:,1);
lon  = input.phase.state(:,2);
lat  = input.phase.state(:,3);
v    = input.phase.state(:,4);
fpa  = input.phase.state(:,5);
azi  = input.phase.state(:,6);
aoa  = input.phase.control(:,1);
bank = input.phase.control(:,2);

% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%
cd0      = input.auxdata.cd(1);
cd1      = input.auxdata.cd(2);
cd2      = input.auxdata.cd(3);
cl0      = input.auxdata.cl(1);
cl1      = input.auxdata.cl(2);
mu       = input.auxdata.mu;
rho0     = input.auxdata.rho0;
H        = input.auxdata.H;
S        = input.auxdata.S;
mass     = input.auxdata.mass;
altitude = rad - input.auxdata.Re;
CD       = cd0+cd1*aoa+cd2*aoa.^2;
rho      = rho0*exp(-altitude/H);
CL       = cl0+cl1*aoa;
q        = 0.5*rho.*v.^2;
D        = q.*S.*CD./mass;
L        = q.*S.*CL./mass;
gravity  = mu./rad.^2;

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%
raddot = v.*sin(fpa);
londot = v.*cos(fpa).*sin(azi)./(rad.*cos(lat));
latdot = v.*cos(fpa).*cos(azi)./rad;
vdot   = -D-gravity.*sin(fpa);
fpadot = (L.*cos(bank)-cos(fpa).*(gravity-v.^2./rad))./v;
azidot = (L.*sin(bank)./cos(fpa)+v.^2.*cos(fpa).*sin(azi).*tan(lat)./rad)./v;

phaseout.dynamics  = [raddot, londot, latdot, vdot, fpadot, azidot];