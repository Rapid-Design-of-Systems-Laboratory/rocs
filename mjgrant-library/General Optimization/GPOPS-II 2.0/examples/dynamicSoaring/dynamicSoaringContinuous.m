%--------------------------------------------%
% BEGIN: function dynamicSoaringContinuous.m %
%--------------------------------------------%
function phaseout = dynamicSoaringContinuous(input)

t             = input.phase(1).time;
s             = input.phase(1).state;
u             = input.phase(1).control;
p             = input.phase(1).parameter;
x             = s(:,1);
y             = s(:,2);
z             = s(:,3);
v             = s(:,4);
gamma         = s(:,5);
psi           = s(:,6);
CL            = u(:,1);
phi           = u(:,2);
beta          = p(:,1);
singamma      = sin(gamma);
cosgamma      = cos(gamma);
sinpsi        = sin(psi);
cospsi        = cos(psi);
sinphi        = sin(phi);
cosphi        = cos(phi);
rho           = input.auxdata.rho;
S             = input.auxdata.S;
CD0           = input.auxdata.CD0;
K             = input.auxdata.K;
g             = input.auxdata.g;
m             = input.auxdata.m;
W0            = input.auxdata.W0;
wx            = (beta.*z+W0);
DWxDt         = beta.*v.*singamma;
vcosgamma     = v.*cosgamma;
DWxDtsinpsi   = DWxDt.*sinpsi;

xdot          = vcosgamma.*sinpsi+wx;
ydot          = vcosgamma.*cospsi;
zdot          = v.*singamma;
term1         = rho*S/2/m;
term2         = 1;
term3         = g*term2;
CLsq          = CL.^2;
vsq           = v.^2;
vdot          = -term1*(CD0+K*CLsq).*vsq-(term3)*singamma-term2*DWxDtsinpsi.*cosgamma;
gammadot      = term1*CL.*v.*cosphi-(term3)*cosgamma./v+term2*DWxDtsinpsi.*singamma./v;
psidot        = (term1*CL.*v.*sinphi-term2*DWxDt.*cospsi./v)./cosgamma;
ngconstant    = (0.5*rho*S/m/g);
ng            = ngconstant.*CL.*v.^2;

phaseout.dynamics  = [xdot, ydot, zdot, vdot, gammadot, psidot];
phaseout.path = ng;

%------------------------------------------%
% END: function dynamicSoaringContinuous.m %
%------------------------------------------%
