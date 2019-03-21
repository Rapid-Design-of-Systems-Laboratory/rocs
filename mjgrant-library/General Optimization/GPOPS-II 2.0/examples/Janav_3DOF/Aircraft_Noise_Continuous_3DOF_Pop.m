function output = Aircraft_Noise_Continuous_3DOF_Pop(input)
%----------------------------------------%
% Begin File:  Aircraft Noise Continuous %
%----------------------------------------%

% Constants

mass = input.auxdata.mass;
g0 = input.auxdata.g0;
C3 = input.auxdata.C3;
C1 = input.auxdata.C1;
C2 = input.auxdata.C2;

% States and time
t = input.phase.time;
state = input.phase.state;
x = state(:,1);
y = state(:,2);
z = state(:,3);
v = state(:,4);
psi = state(:,5);
gam = state(:,6);

% Controls
bank = input.phase.control(:,1);
aoa  = input.phase.control(:,2);
T = input.phase.control(:,3);

% Important parameters
L = mass.*g0; % L = W
D = 0.226.*v.^2 + 5.2e6./(v.^2);
gravity = g0;
sgam = sin(gam);
cgam = cos(gam);

% Approximate population
% keyboard
% for count = 1:1:length(x)
% pop(count) = population(x(count)/1000,y(count)/1000);
% end
% pop = -0.093425*x.^3 + 0.017248*x.^2.*y +...
%     0.040138*x.^2 - 0.0035582*x.*y.^2 + 0.076212*x.*y +...
%     1.3537*x - 0.0069758*x.^3 - 0.059675*y.^2 + 0.062387*y + 2.481;

% Below code works
% pop = 1;
%

% integrand  = log(C3)+5.2*log(T)+log(cos(gam))-(log(v)+2.5*log(h+50));
% x1 = x/1000;
% y1 = y/1000;
pop = exp(5500-x);
% pop = 2*(5500 -x);
% pop = -0.093425*x1.^3 + 0.017248.*x1.^2.*y1 + 0.74082.*x1.^2 - 0.0035582.*x1.*y1.^2 - 0.010029.*x1.*y1 - 0.59874.*x1 - 0.0069758.*y1.^3 - 0.05078.*y1.^2 - 0.020341.*y1 + 0.8075;
integrand = pop.*C3.*T.^(5.2).*cos(gam)./(v.*(z+50).^(2.5));
% integrand = C3*T.^2./(v.*(h+50));

% Equations of motion
xdot = v.*cgam.*cos(psi);
ydot = v.*cgam.*sin(psi);
zdot   = v.*sgam;
vdot = (T.*cos(aoa)-D)./mass-gravity.*sgam;
psidot = (T.*sin(aoa)+L).*sin(bank)./(mass.*v.*cgam);
gamdot   = ((T.*sin(aoa)+L).*cos(bank)./mass-cgam.*gravity)./v;
output.dynamics  = [xdot, ydot, zdot ,vdot, psidot, gamdot];
output.integrand = integrand;
%------------------------------------%
% End File:  goddardRocketContinuous %
%------------------------------------%
