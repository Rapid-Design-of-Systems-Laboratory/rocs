%-------------------------------------------%
% Begin Function:  tuberculosisContinuous.m %
%-------------------------------------------%
function phaseout = tuberculosisContinuous(input)

beta1 = input.auxdata.beta1;
beta2 = input.auxdata.beta2;
mu = input.auxdata.mu;
d1 = input.auxdata.d1;
d2 = input.auxdata.d2;
k1 = input.auxdata.k1;
k2 = input.auxdata.k2;
r1 = input.auxdata.r1;
r2 = input.auxdata.r2;
p = input.auxdata.p;
q = input.auxdata.q;
Npop = input.auxdata.Npop;
betas = input.auxdata.betas;
B1 = input.auxdata.B1;
B2 = input.auxdata.B2;
lam = input.auxdata.lam;

S = input.phase.state(:,1);
T  = input.phase.state(:,2);
L1 = input.phase.state(:,3);
L2 = input.phase.state(:,4);
I1 = input.phase.state(:,5);
I2 = input.phase.state(:,6);

u1 = input.phase.control(:,1);
u2 = input.phase.control(:,2);

dS = lam-(beta1.*S.*I1+betas.*S.*I2)./Npop-mu.*S;
dT = u1.*r1.*L1-mu.*T+(1-(1-u2).*(p+q)).*r2.*I1-(beta2.*T.*I1+betas.*T.*I2)./Npop;
dL1 = (beta1.*S.*I1+beta2.*T.*I1-betas.*L1.*I2)./Npop...
      -(mu+k1).*L1-u1.*r1.*L1+(1-u2).*p.*r2.*I1;
dL2 = (1-u2).*q.*r2.*I1-(mu+k2).*L2+betas.*(S+L1+T).*I2./Npop;
dI1 = k1.*L1-(mu+d1).*I1-r2.*I1;
dI2 = k2.*L2-(mu+d2).*I2;

phaseout.dynamics  = [dS, dT, dL1, dL2, dI1, dI2];
phaseout.path = S + T + L1 + L2 + I1 + I2 - Npop;
phaseout.integrand = L2 + I2 + B1./2.*u1.^2 + B2./2.*u2.^2;

%-------------------------------------------%
% End Function:  tuberculosisContinuous.m   %
%-------------------------------------------%
