function XDot = dynamics(t,x,auxdata,guess)
% ---------------------------------------------------------%
% This function is used to propagate the dynamics in order %
% to obtain an initial guess of the solution for GPOPS-II  %
% ---------------------------------------------------------%

T = auxdata.T;
Isp = auxdata.Isp;
mu = auxdata.mu;
gs = auxdata.gs;
Re = auxdata.Re;
J2 = auxdata.J2;
J3 = auxdata.J3;
J4 = auxdata.J4;

tau = guess.tau;

p = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6);
w = x(7);

% ----------------------------------------------------------------------- %
% Gravitational disturbing acceleration
% ----------------------------------------------------------------------- %
q = 1+f*cos(L)+g*sin(L);
r = p/q;
alpha2 = h*h-k*k;
chi = sqrt(h*h+k*k);
s2 = 1+chi*chi;

rX = (r/s2)*(cos(L)+alpha2*cos(L)+2*h*k*sin(L));
rY = (r/s2)*(sin(L)-alpha2*sin(L)+2*h*k*cos(L));
rZ = (2*r/s2)*(h*sin(L)-k*cos(L));
rVec = [rX rY rZ];
rMag = sqrt(rX^2+rY^2+rZ^2);
rXZMag = sqrt(rX^2+rZ^2);

vX = -(1/s2)*sqrt(mu/p)*(sin(L)+alpha2*sin(L)-2*h*k*cos(L)+g-2*f*h*k+alpha2*g);
vY = -(1/s2)*sqrt(mu/p)*(-cos(L)+alpha2*cos(L)+2*h*k*sin(L)-f+2*g*h*k+alpha2*f);
vZ = (2/s2)*sqrt(mu/p)*(h*cos(L)+k*sin(L)+f*h+g*k);
vVec = [vX vY vZ];
vMag = sqrt(vX^2+vY^2+vZ^2);

rCrossv = cross(rVec,vVec,2);
rCrossvMag = sqrt(rCrossv(1)^2+rCrossv(2)^2+rCrossv(3)^2);
rCrossvCrossr = cross(rCrossv,rVec,2);

ir1 = rVec(1)/rMag;
ir2 = rVec(2)/rMag;
ir3 = rVec(3)/rMag;
ir = [ir1 ir2 ir3];

it1 = rCrossvCrossr(1)/(rCrossvMag*rMag);
it2 = rCrossvCrossr(2)/(rCrossvMag*rMag);
it3 = rCrossvCrossr(3)/(rCrossvMag*rMag);
it = [it1 it2 it3];

ih1 = rCrossv(1)/rCrossvMag;
ih2 = rCrossv(2)/rCrossvMag;
ih3 = rCrossv(3)/rCrossvMag;
ih = [ih1 ih2 ih3];

enir = ir3;

enirir1 = enir.*ir1;
enirir2 = enir.*ir2;
enirir3 = enir.*ir3;

enenirir1 = 0-enirir1;
enenirir2 = 0-enirir2;
enenirir3 = 1-enirir3;
enenirirMag = sqrt(enenirir1^2+enenirir2^2+enenirir3^2);

in1 = enenirir1/enenirirMag;
in2 = enenirir2/enenirirMag;
in3 = enenirir3/enenirirMag;

% Geocentric latitude angle
sinphi = rZ/rXZMag;
cosphi = sqrt(1-sinphi^2);

% Legendre polynomials
P2 = (3*sinphi^2-2)/2;
P3 = (5*sinphi^3-3*sinphi)/2;
P4 = (35*sinphi^4-30*sinphi^2+3)/8;
dP2 = 3*sinphi;
dP3 = (15*sinphi-3)/2;
dP4 = (140*sinphi^3-60*sinphi)/8;

% Oblate earth perturbations
sumn = (Re/r)^2*dP2*J2+(Re/r)^3*dP3*J3+(Re/r)^4*dP4*J4;
sumr = (2+1)*(Re/r)^2*P2*J2+(3+1)*(Re/r)^3*P3*J3+(4+1)*(Re/r)^4*P4*J4;
deltagn = -(mu*cosphi/(r^2))*sumn;
deltagr = -(mu/(r^2))*sumr;

deltagnin1 = deltagn*in1;
deltagnin2 = deltagn*in2;
deltagnin3 = deltagn*in3;

deltagrir1 = deltagr*ir1;
deltagrir2 = deltagr*ir2;
deltagrir3 = deltagr*ir3;

deltag1 = deltagnin1 - deltagrir1;
deltag2 = deltagnin2 - deltagrir2;
deltag3 = deltagnin3 - deltagrir3;

Deltag1 = ir(1)*deltag1+ir(2)*deltag2+ir(3)*deltag3;
Deltag2 = it(1)*deltag1+it(2)*deltag2+it(3)*deltag3;
Deltag3 = ih(1)*deltag1+ih(2)*deltag2+ih(3)*deltag3;

Qr = [ir' it' ih'];

vUnit = (vVec/vMag)';

u = Qr'*vUnit;
ur = u(1);
ut = u(2);
uh = u(3);

% ----------------------------------------------------------------------- %
% Thrust acceleration
% ----------------------------------------------------------------------- %
DeltaT1 = ((gs*T*(1+0.01*tau))/w)*ur;
DeltaT2 = ((gs*T*(1+0.01*tau))/w)*ut;
DeltaT3 = ((gs*T*(1+0.01*tau))/w)*uh;

% ----------------------------------------------------------------------- %
% Total acceleration
% ----------------------------------------------------------------------- %
Delta1 = Deltag1+DeltaT1;
Delta2 = Deltag2+DeltaT2;
Delta3 = Deltag3+DeltaT3;

q = 1+f*cos(L)+g*sin(L);
s2 = 1+h^2+k^2;

dp = (2*p/q)*sqrt(p/mu)*Delta2;

df =  sqrt(p/mu)*sin(L)*Delta1 ...
     +sqrt(p/mu)*(1/q)*((q+1)*cos(L)+f)*Delta2 ...
     -sqrt(p/mu)*(g/q)*(h*sin(L)-k*cos(L))*Delta3;
 
dg = -sqrt(p/mu)*cos(L)*Delta1 ...
     +sqrt(p/mu)*(1/q)*((q+1)*sin(L)+g)*Delta2 ...
     +sqrt(p/mu)*(f/q)*(h*sin(L)-k*cos(L))*Delta3;
 
dh = sqrt(p/mu)*(s2*cos(L)/(2*q))*Delta3;

dk = sqrt(p/mu)*(s2*sin(L)/(2*q))*Delta3;

dL =  sqrt(p/mu)*(1/q)*(h*sin(L)-k*cos(L))*Delta3...
     +sqrt(mu*p)*((q/p)^2);
 
dw = -(T*(1+0.01*tau)/Isp);

XDot = [dp; df; dg; dh; dk; dL; dw];

end
