function phaseout = seywaldMinimumTimeToClimbContinuous(input)

g    = input.auxdata.g;
m    = input.auxdata.m;
S    = input.auxdata.S;
rho0 = input.auxdata.rho0;
a0   = input.auxdata.a0;

t   = input.phase.time;
h   = input.phase.state(:,1);
v   = input.phase.state(:,2);
fpa = input.phase.state(:,3);
u   = input.phase.control;


  hbar = h./1000;
  HH   = zeros(size(t,1),5);
  HH(:,1) = ones(size(t));
  for ii=2:6
    HH(:,ii) = HH(:,ii-1).*hbar; 
  end;
  z = HH(:,2:5)*input.auxdata.z;
  r = input.auxdata.r0*exp(-z);
  y = HH(:,1:2)*input.auxdata.y + r;
  rho = rho0*exp(y);
  theta = HH(:,1:4)*input.auxdata.theta;
  a     = a0*sqrt(theta);
  
  M = v./a;
  
  q = rho.*v.*v.*input.auxdata.S/2;
  L = m*g*u;
  
  MM   = zeros(size(t,1),6);
  MM(:,1) = ones(size(t));
  for ii=2:6
    MM(:,ii) = MM(:,ii-1).*M; 
  end;
  numeratorCD0   = MM(:,1:6)*input.auxdata.a;
  denominatorCD0 = MM(:,1:6)*input.auxdata.b;
  CD0            = numeratorCD0./denominatorCD0;
  numeratorK     = MM(:,1:6)*input.auxdata.c;
  denominatorK   = MM(:,1:6)*input.auxdata.d;
  K              = numeratorK./denominatorK;
  D = q.*(CD0+K.*((m.^2).*(g.^2)./(q.^2)).*(u.^2));
  
  E = MM(:,1:6)*input.auxdata.f;
  % T  = sum(HH(:,1:6).*E,2).*g/2.2;
  T  = sum(HH(:,1:6).*E,2);
  
  hdot = v.*sin(fpa);
  vdot = (T-D)./m-g.*sin(fpa);
  fpadot = g.*(u-cos(fpa))./v;
  

phaseout(1).dynamics = [hdot, vdot, fpadot];

%---------------------------------%
% END: function minimumClimbDae.m %
%---------------------------------%
