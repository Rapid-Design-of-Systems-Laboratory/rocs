function phaseout = spaceStationContinuous(input)

time = input.phase.time;
omega = input.phase.state(:,1:3);
r = input.phase.state(:,4:6);
h = input.phase.state(:,7:9);
u=input.phase.control;
N = length(time);

J=input.auxdata.J;
omegaorb = input.auxdata.omegaorb;

hmag = sqrt(dot(h,h,2));
hmagsq = dot(h,h,2);

%%|Calculating Omegadot|

a1=(J(1,1)*omega(:,1))+J(1,2)*omega(:,2)+J(1,3)*omega(:,3)+h(:,1); %column 1 of A
a2=(J(2,1)*omega(:,1))+J(2,2)*omega(:,2)+J(2,3)*omega(:,3)+h(:,2); %column 2 of A
a3=(J(3,1)*omega(:,1))+J(3,2)*omega(:,2)+J(3,3)*omega(:,3)+h(:,3); %column 3 of A

A=[a1,a2,a3];
B=cross(omega,A,2);

%|Calculating C|

D=2./(1+dot(r,r,2));

%r cross-Identity
e1=[-ones(N,1),r(:,3),-r(:,2)];
e2=[-r(:,3),-ones(N,1),r(:,1)];
e3=[r(:,2),-r(:,1),-ones(N,1)];

%r cross E
f1=cross(r,e1,2);
f2=cross(r,e2,2);
f3=cross(r,e3,2);

%E+Identity
C1=[f1(:,1).*D+ones(N,1),f1(:,2).*D,f1(:,3).*D];
C2=[f2(:,1).*D,f2(:,2).*D+ones(N,1),f2(:,3).*D];
C3=[f3(:,1).*D,f3(:,2).*D,f3(:,3).*D+ones(N,1)];

I=[C3*J]; %------------------------------------------------------

C3crossJC3=cross(C3,I,2);
threeomegaorbsq=3*(omegaorb^2);

taugg=threeomegaorbsq*C3crossJC3;

%term multiplied by J inverse to get omegadot
K=taugg-B-u;

Jinv = input.auxdata.Jinv;

omegadot=K*Jinv;   %omegadot   ----------------------------------

%%|Calculating rdot|

omegaor=-omegaorb.*C2;

omegaminusomegaor = omega-omegaor;

%r times r transpose + Identity + rcross
l1=[r(:,1).^2+1,r(:,2).*r(:,1)+r(:,3),r(:,3).*r(:,1)-r(:,2)];    %--------
l2=[r(:,1).*r(:,2)-r(:,3),r(:,2).^2+1,r(:,3).*r(:,2)+r(:,1)];    %--------
l3=[r(:,1).*r(:,3)+r(:,2),r(:,2).*r(:,3)-r(:,1),r(:,3).^2+1];    %--------

rdot1=0.5*[l1(:,1).*omegaminusomegaor(:,1)+l2(:,1).*omegaminusomegaor(:,2)+l3(:,1).*omegaminusomegaor(:,3)];  %rows
rdot2=0.5*[l1(:,2).*omegaminusomegaor(:,1)+l2(:,2).*omegaminusomegaor(:,2)+l3(:,2).*omegaminusomegaor(:,3)];
rdot3=0.5*[l1(:,3).*omegaminusomegaor(:,1)+l2(:,3).*omegaminusomegaor(:,2)+l3(:,3).*omegaminusomegaor(:,3)];

rdot = [rdot1,rdot2,rdot3];   %----------------------------------------- 

hdot=u;  %hdot

phaseout.dynamics  = [omegadot, rdot, hdot];
phaseout.path      = hmagsq;
phaseout.integrand = 1e-6*dot(u,u,2);
