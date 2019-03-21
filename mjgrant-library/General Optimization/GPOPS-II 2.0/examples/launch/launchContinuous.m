%-------------------------------------------------------------------------%
%------------------- BEGIN Function launchContinuous.m -------------------%
%-------------------------------------------------------------------------%
function phaseout = launchContinuous(input)

%---------------------%
% Dynamics in Phase 1 %
%---------------------%
t1 = input.phase(1).time;
x1 = input.phase(1).state;
u1 = input.phase(1).control;
r1 = x1(:,1:3);
v1 = x1(:,4:6);
m1 = x1(:,7);

rad1 = sqrt(sum(r1.*r1,2));
omegaMatrix = input.auxdata.omegaMatrix;
omegacrossr = r1*omegaMatrix.';
vrel1 = v1-omegacrossr;
speedrel1 = sqrt(sum(vrel1.*vrel1,2));
h1 = rad1-input.auxdata.Re;
rho1 = input.auxdata.rho0*exp(-h1/input.auxdata.H);
bc1  = (rho1./(2*m1)).*input.auxdata.sa*input.auxdata.cd;
bcspeed1 = bc1.*speedrel1;
bcspeedmat1 = repmat(bcspeed1,1,3);
Drag1 = -bcspeedmat1.*vrel1;
muoverradcubed1 = input.auxdata.mu./rad1.^3;
muoverradcubedmat1 = [muoverradcubed1 muoverradcubed1 muoverradcubed1];
grav1 = -muoverradcubedmat1.*r1;

TSrb1   = 6*input.auxdata.thrustSrb*ones(size(t1));
TFirst1 = input.auxdata.thrustFirst*ones(size(t1));
TTot1   = TSrb1+TFirst1;
m1dot1  = -TSrb1./(input.auxdata.g0*input.auxdata.ispSrb);
m2dot1  = -TFirst1./(input.auxdata.g0*input.auxdata.ispFirst);
mdot1   = m1dot1+m2dot1;
q1 = 1/2*rho1.*speedrel1.^2;
path1 = [sum(u1.*u1,2)];
Toverm1 = TTot1./m1;
Tovermmat1 = [Toverm1 Toverm1 Toverm1];
thrust1 = Tovermmat1.*u1;
rdot1 = v1;
vdot1 = thrust1+Drag1+grav1;
phaseout(1).dynamics = [rdot1 vdot1 mdot1];
phaseout(1).path = path1;

%---------------------%
% Dynamics in Phase 2 %
%---------------------%
t2 = input.phase(2).time;
x2 = input.phase(2).state;
u2 = input.phase(2).control;
r2 = x2(:,1:3);
v2 = x2(:,4:6);
m2 = x2(:,7);

rad2 = sqrt(sum(r2.*r2,2));
omegaMatrix = input.auxdata.omegaMatrix;
omegacrossr = r2*omegaMatrix.';
vrel2 = v2-omegacrossr;
speedrel2 = sqrt(sum(vrel2.*vrel2,2));
h2 = rad2-input.auxdata.Re;
rho2 = input.auxdata.rho0*exp(-h2/input.auxdata.H);
bc2  = (rho2./(2*m2)).*input.auxdata.sa*input.auxdata.cd;
bcspeed2 = bc2.*speedrel2;
bcspeedmat2 = repmat(bcspeed2,1,3);
Drag2 = -bcspeedmat2.*vrel2;
muoverradcubed2 = input.auxdata.mu./rad2.^3;
muoverradcubedmat2 = [muoverradcubed2 muoverradcubed2 muoverradcubed2];
grav2 = -muoverradcubedmat2.*r2;
TSrb2 = 3*input.auxdata.thrustSrb*ones(size(t2));
TFirst2 = input.auxdata.thrustFirst*ones(size(t2));
TTot2 = TSrb2+TFirst2;
m1dot2 = -TSrb2./(input.auxdata.g0*input.auxdata.ispSrb);
m2dot2 = -TFirst2./(input.auxdata.g0*input.auxdata.ispFirst);
mdot2 = m1dot2+m2dot2;    
path2 = [sum(u2.*u2,2)];
Toverm2 = TTot2./m2;
Tovermmat2 = [Toverm2 Toverm2 Toverm2];
thrust2 = Tovermmat2.*u2;
rdot2 = v2;
vdot2 = thrust2+Drag2+grav2;
phaseout(2).dynamics = [rdot2 vdot2 mdot2];
phaseout(2).path = path2;

%---------------------%
% Dynamics in Phase 3 %
%---------------------%
t3 = input.phase(3).time;
x3 = input.phase(3).state;
u3 = input.phase(3).control;
r3 = x3(:,1:3);
v3 = x3(:,4:6);
m3 = x3(:,7);

rad3 = sqrt(sum(r3.*r3,2));
omegaMatrix = input.auxdata.omegaMatrix;
omegacrossr = r3*omegaMatrix.';
vrel3 = v3-omegacrossr;
speedrel3 = sqrt(sum(vrel3.*vrel3,2));
h3 = rad3-input.auxdata.Re;
rho3 = input.auxdata.rho0*exp(-h3/input.auxdata.H);
bc3  = (rho3./(2*m3)).*input.auxdata.sa*input.auxdata.cd;
bcspeed3 = bc3.*speedrel3;
bcspeedmat3 = repmat(bcspeed3,1,3);
Drag3 = -bcspeedmat3.*vrel3;
muoverradcubed3 = input.auxdata.mu./rad3.^3;
muoverradcubedmat3 = [muoverradcubed3 muoverradcubed3 muoverradcubed3];
grav3 = -muoverradcubedmat3.*r3;
TTot3 = input.auxdata.thrustFirst*ones(size(t3));
mdot3 = -TTot3./(input.auxdata.g0*input.auxdata.ispFirst);
path3 = [sum(u3.*u3,2)];
Toverm3 = TTot3./m3;
Tovermmat3 = [Toverm3 Toverm3 Toverm3];
thrust3 = Tovermmat3.*u3;
rdot3 = v3;
vdot3 = thrust3+Drag3+grav3;
phaseout(3).dynamics = [rdot3 vdot3 mdot3];
phaseout(3).path = path3;

%---------------------%
% Dynamics in Phase 4 %
%---------------------%
t4 = input.phase(4).time;
x4 = input.phase(4).state;
u4 = input.phase(4).control;
r4 = x4(:,1:3);
v4 = x4(:,4:6);
m4 = x4(:,7);
rad4 = sqrt(sum(r4.*r4,2));
omegacrossr = r4*input.auxdata.omegaMatrix.';
vrel4 = v4-omegacrossr;
speedrel4 = sqrt(sum(vrel4.*vrel4,2));
h4 = rad4-input.auxdata.Re;
rho4 = input.auxdata.rho0*exp(-h4/input.auxdata.H);
bc4  = (rho4./(2*m4)).*input.auxdata.sa*input.auxdata.cd;
bcspeed4 = bc4.*speedrel4;
bcspeedmat4 = repmat(bcspeed4,1,3);
Drag4 = -bcspeedmat4.*vrel4;
muoverradcubed4 = input.auxdata.mu./rad4.^3;
muoverradcubedmat4 = [muoverradcubed4 muoverradcubed4 muoverradcubed4];
grav4 = -muoverradcubedmat4.*r4;
TTot4 = input.auxdata.thrustSecond*ones(size(t4));
mdot4 = -TTot4/(input.auxdata.g0*input.auxdata.ispSecond);
path4 = [sum(u4.*u4,2)];
Toverm4 = TTot4./m4;
Tovermmat4 = [Toverm4 Toverm4 Toverm4];
thrust4 = Tovermmat4.*u4;
rdot4 = v4;
vdot4 = thrust4+Drag4+grav4;
phaseout(4).dynamics = [rdot4 vdot4 mdot4];
phaseout(4).path = path4;

%-------------------------------------------------------------------------%
%--------------------- END Function launchContinuous.m -------------------%
%-------------------------------------------------------------------------%
