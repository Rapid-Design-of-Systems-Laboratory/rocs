function initialguess = lowThrustPropagate(auxdata,initial,guess)

% ----------------------------------------------------------------------- %
% Constants and initial conditions
% ----------------------------------------------------------------------- %
T = auxdata.T;
Isp = auxdata.Isp;
mu = auxdata.mu;
gs = auxdata.gs;
Re = auxdata.Re;
J2 = auxdata.J2;
J3 = auxdata.J3;
J4 = auxdata.J4;

p0 = initial.p0; % [ft]
f0 = initial.f0;
g0 = initial.g0;
h0 = initial.h0;
k0 = initial.k0;
L0 = initial.L0;
w0 = initial.w0;
t0 = initial.t0;

tau = guess.tau;
tf  = guess.tf;

% ----------------------------------------------------------------------- %
% Integrate ordinary differential equations using ode45
% ----------------------------------------------------------------------- %
TSPAN = [t0,tf];
X0 = [p0; f0; g0; h0; k0; L0; w0];
ODEOPTS = odeset('RelTol',1e-6,'AbsTol',1e-6);

tic;
[TOUT,XOUT] = ode113(@lowThrustDynamics,TSPAN,X0,ODEOPTS,auxdata,guess);
propCPU = toc;

% ----------------------------------------------------------------------- %
% Extract solution
% ----------------------------------------------------------------------- %
p = XOUT(:,1);
f = XOUT(:,2);
g = XOUT(:,3);
h = XOUT(:,4);
k = XOUT(:,5);
L = XOUT(:,6);
w = XOUT(:,7);

% Convert to Cartesian
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;

rX = (r./s2).*(cos(L)+alpha2.*cos(L)+2*h.*k.*sin(L));
rY = (r./s2).*(sin(L)-alpha2.*sin(L)+2*h.*k.*cos(L));
rZ = (2*r./s2).*(h.*sin(L)-k.*cos(L));
rVec = [rX rY rZ];
rMag = sqrt(rX.^2+rY.^2+rZ.^2);

vX = -(1./s2).*sqrt(mu./p).*(sin(L)+alpha2.*sin(L)-2*h.*k.*cos(L)+g-2*f.*h.*k+alpha2.*g);
vY = -(1./s2).*sqrt(mu./p).*(-cos(L)+alpha2.*cos(L)+2*h.*k.*sin(L)-f+2*g.*h.*k+alpha2.*f);
vZ = (2./s2).*sqrt(mu./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);
vVec = [vX vY vZ];
vMag = sqrt(vX.^2+vY.^2+vZ.^2);

% Calculate control components
rCrossv = cross(rVec,vVec,2);
rCrossvMag = sqrt(rCrossv(:,1).^2+rCrossv(:,2).^2+rCrossv(:,3).^2);
rCrossvCrossr = cross(rCrossv,rVec,2);

ir1 = rVec(:,1)./rMag;
ir2 = rVec(:,2)./rMag;
ir3 = rVec(:,3)./rMag;
ir = [ir1 ir2 ir3];

it1 = rCrossvCrossr(:,1)./(rCrossvMag.*rMag);
it2 = rCrossvCrossr(:,2)./(rCrossvMag.*rMag);
it3 = rCrossvCrossr(:,3)./(rCrossvMag.*rMag);
it = [it1 it2 it3];

ih1 = rCrossv(:,1)./rCrossvMag;
ih2 = rCrossv(:,2)./rCrossvMag;
ih3 = rCrossv(:,3)./rCrossvMag;
ih = [ih1 ih2 ih3];

vUnit1 = vVec(:,1)./vMag;
vUnit2 = vVec(:,2)./vMag;
vUnit3 = vVec(:,3)./vMag;

ur = ir(:,1).*vUnit1+ir(:,2).*vUnit2+ir(:,3).*vUnit3;
ut = it(:,1).*vUnit1+it(:,2).*vUnit2+it(:,3).*vUnit3;
uh = ih(:,1).*vUnit1+ih(:,2).*vUnit2+ih(:,3).*vUnit3;

% Convert MEE to COE
[a,ecc,inc,Ome,ome,nu] = lowThrustMee2Coe(p,f,g,h,k,L);

% Declare initial guess for GPOPS-II
pguess = p;
fguess = f;
gguess = g;
hguess = h;
kguess = k;
Lguess = L;
wguess = w;
urguess = ur;
uhguess = uh;
utguess = ut;
tauguess = tau;
tguess = TOUT;

% ----------------------------------------------------------------------- %
% Save data to *.mat file
% ----------------------------------------------------------------------- %
% save('propdata.mat','pguess','fguess','gguess','hguess','kguess','Lguess','wguess','urguess','uhguess','utguess','tauguess','tguess','propCPU');

% Determine initial mesh for GPOPS2
M = ceil(L(end)/(2*pi));
n = 1:1:M;
Lmesh = L(1)+((n-1)./(M-1)).*(L(end)-L(1));
tmesh = spline(Lguess,TOUT,Lmesh);
fraction = zeros(1,length(tmesh)-1);
for i = 1:length(tmesh)-1   
    fraction(1,i) = tmesh(1,i+1)-tmesh(1,i);
end
fraction = fraction./sum(fraction);

% ----------------------------------------------------------------------- %
% Save initial mesh setup to *.mat file
% ----------------------------------------------------------------------- %
initialguess.phase.time      = tguess;
initialguess.phase.state     = [pguess,fguess,gguess,hguess,kguess,Lguess,wguess];
initialguess.phase.control   = [urguess utguess uhguess];
initialguess.parameter       = guess.tau;
initialguess.fraction        = fraction;
