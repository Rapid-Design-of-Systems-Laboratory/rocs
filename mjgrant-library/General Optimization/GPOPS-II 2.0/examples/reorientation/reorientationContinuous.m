function phaseout = reorientationContinuous(input)

% input
% input.phase(phasenumber).state
% input.phase(phasenumber).control
% input.phase(phasenumber).time
% input.phase(phasenumber).parameter
%
% input.auxdata = auxiliary information
%
% output
% phaseout(phasenumber).dynamics
% phaseout(phasenumber).path
% phaseout(phasenumber).integrand

% Time
t      = input.phase.time;

% States
q1     = input.phase.state(:,1);
q2     = input.phase.state(:,2);
q3     = input.phase.state(:,3);
w1     = input.phase.state(:,4);
w2     = input.phase.state(:,5);
w3     = input.phase.state(:,6);

% Control
q4 = input.phase.control(:,1);
u1 = input.phase.control(:,2);
u2 = input.phase.control(:,3);
u3 = input.phase.control(:,4);

% Constants being used
   Ix = 5621; Iy = 4547; Iz = 2364;
% Differential Equations for the states
  dq1 = (1/2).*( w1.*q4 - w2.*q3 + w3.*q2);              
  dq2 = (1/2).*( w1.*q3 + w2.*q4 - w3.*q1);              
  dq3 = (1/2).*(-w1.*q2 + w2.*q1 + w3.*q4);             
  dw1 = (u1./Ix)-((Iz-Iy)/Ix).*w2.*w3;            
  dw2 = (u2./Iy)-((Ix-Iz)/Iy).*w1.*w3;           
  dw3 = (u3./Iz)-((Iy-Ix)/Iz).*w1.*w2;    

phaseout.dynamics  = [dq1, dq2, dq3, dw1, dw2, dw3]; 
phaseout.path = q1.^2 + q2.^2 + q3.^2 + q4.^2;
