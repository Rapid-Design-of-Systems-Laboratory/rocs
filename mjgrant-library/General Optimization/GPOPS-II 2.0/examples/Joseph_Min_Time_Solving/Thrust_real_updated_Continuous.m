function output = Thrust_real_updated_Continuous(input)
%--------------------------------------%
% Begin File:  HMME_phases_Continuous %
%--------------------------------------%

% auxdata = input.auxdata;
cd0      = input.auxdata.cd(1);
cd1      = input.auxdata.cd(2);
cd2      = input.auxdata.cd(3);
cl0      = input.auxdata.cl(1);
cl1      = input.auxdata.cl(2);

% keyboard 

mu  = input.auxdata.mu;
rho0 = input.auxdata.rho0;
H = input.auxdata.H;
S = input.auxdata.S;
% Isp = input.auxdata.Isp; 
g0 = input.auxdata.g0;

% New parameters
y  = input.auxdata.y;
R  = input.auxdata.R;
T0  = input.auxdata.T0;
T4  = input.auxdata.T4;
cp  = input.auxdata.cp;
hpr  = input.auxdata.hpr;
Mc = input.auxdata.Mc;
M1 = input.auxdata.M1; 
M2 = input.auxdata.M2;  
q1 = input.auxdata.q1;  
q2 = input.auxdata.q2; 
qscale = input.auxdata.qscale;
T1 = input.auxdata.T1; 
t = input.phase.time;
x = input.phase.state;
h = x(:,1);
thetta = x(:,2);
v = x(:,3);
gam = x(:,4);
mass = x(:,5);
aoa = x(:,6);
alfadot  = input.phase.control(:,1);
A = input.phase.control(:,2); % This is new control
Re = input.auxdata.Re;

% Calculate stuff for phase 1
rad = h + Re;
rho = rho0*exp(-h/H);
CD       = cd0+cd1*aoa+cd2*aoa.^2;
CL       = cl0+cl1*aoa;
dynamic_pressure = 0.5*rho.*v.^2;
q = 0.5*rho.*v.^2;
D = dynamic_pressure.*S.*CD;
L = dynamic_pressure.*S.*CL;
gravity = mu./rad.^2;
bank = 0;
sgam = sin(gam);
cgam = cos(gam);
cbank = cos(bank);

% New equations needed
va = sqrt(y*R*T0);
M = v./va;
rho = rho0*exp(-h/H); % Exponential Atmospheric Density [kg/m^3]
madot = rho.*A.*v; % air flow rate
Tau_lam = T4*(1+(y-1)*Mc^2/2)/T0;
Tau_r = (1+(y-1).*M.^2/2);
%Tau_lam = Tau_r +(Tau_lam - Tau_r).*(atan(qscale*(q-q1))-atan(qscale*(q-q2)))/pi;
f = cp*T0/hpr*(Tau_lam-Tau_r); % fuel ratio
ST = v.*(sqrt(Tau_lam./Tau_r)-1);
Isp = ST./(f*g0);
% % keyboard
% for i = 1:1:length(ST)
%     % M_info = M(i);
%     q_info = dynamic_pressure(i);
%     if (q_info<40000 || q_info>150000) % M_info<4 || M_info>9.75 || 
%         T(i) = 0;
%     else
%       T(i) = madot(i)*ST(i);  
%     end
% end
% T = T';

T = madot.*ST;
% T = T.*(atan(1e17*(M-M1))+atan(1e17*(M-M2)))/pi;

% Tmax = input.auxdata.Tmax;
% T = Tmax*A;
% path = dynamic_pressure;
% output.path = path;

hdot   = v.*sgam;
thettadot = v.*cgam./rad;
vdot = (T.*cos(aoa)-D)./mass-gravity.*sgam;
gamdot   = ((T.*sin(aoa)+L.*cbank)./mass-cgam.*(gravity-v.^2./rad))./v;
massdot = -madot.*f;
aoadot = alfadot;
output.dynamics  = [hdot, thettadot, vdot, gamdot, massdot, aoadot];

%------------------------------------%
% End File:  goddardRocketContinuous %
%------------------------------------%
