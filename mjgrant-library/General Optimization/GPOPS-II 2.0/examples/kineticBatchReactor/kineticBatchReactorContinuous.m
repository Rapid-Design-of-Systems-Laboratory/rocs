%--------------------------------------%
% BEGIN: kineticBatchReactorContinuous.m
%--------------------------------------%
function phaseout = kineticBatchReactorContinuous(input)

%---------------------------------------%
% EXTRACT THE DATA COMMON TO ALL PHASES %
%---------------------------------------%
auxdata    = input.auxdata;
beta1      = auxdata.beta1;
betaminus1 = auxdata.betaminus1;
beta2      = auxdata.beta2;
K1         = auxdata.K1;
K2         = auxdata.K2;
K3         = auxdata.K3;

%---------------------------------------------------------------------%
% Get the Differential Equations and Path Constraints in Phase IPHASE %
%---------------------------------------------------------------------%
for iphase=1:3

%---------------------------------------------------------%
% Extract the Physical Data Corresponding to Phase IPHASE %
%---------------------------------------------------------%
k1hat      = auxdata.phase(iphase).k1hat;
kminus1hat = auxdata.phase(iphase).kminus1hat;
k2hat      = auxdata.phase(iphase).k2hat;
a          = auxdata.phase(iphase).a;

%-----------------------------------------------------------------%
% Extract the Time, State, Control and Parameters in Phase IPHASE %
%-----------------------------------------------------------------%
t          = input.phase(iphase).time;
y1         = input.phase(iphase).state(:,1);
y2         = input.phase(iphase).state(:,2);
y3         = input.phase(iphase).state(:,3);
y4         = input.phase(iphase).state(:,4);
y5         = input.phase(iphase).state(:,5);
y6         = input.phase(iphase).state(:,6);
u1         = input.phase(iphase).control(:,1);
u2         = input.phase(iphase).control(:,2);
u3         = input.phase(iphase).control(:,3);
u4         = input.phase(iphase).control(:,4);
u5         = input.phase(iphase).control(:,5);
p          = input.phase(iphase).parameter;
k1         = k1hat*exp(-beta1./u5);
kminus1    = kminus1hat*exp(-betaminus1./u5);
k2         = k2hat*exp(-beta2./u5);
k3         = k1;
kminus3    = 0.5*kminus1;

%----------------------------------------------------%
% Compute the Differential Equations in Phase IPHASE %
%----------------------------------------------------%
y1dot      = -k2.*y2.*u2;
y2dot      = -k1.*y2.*y6+kminus1.*u4+y1dot;
y3dot      = -y1dot+k3.*y4.*y6-kminus3.*u3;
y4dot      = -k3.*y4.*y6+kminus3.*u3;
y5dot      =  k1.*y2.*y6-kminus1.*u4;
y6dot      = -k1.*y2.*y6+kminus1.*u4-k3.*y4.*y6+kminus3.*u3;
path1      =  1*(p-y6+10.^(-u1)-u2-u3-u4);	 	
path2      =  1*(u2-K2*y1./(K2+10.^(-u1)));
path3      =  1*(u3-K3*y3./(K3+10.^(-u1)));	
path4      =  1*(u4-K1*y5./(K1+10.^(-u1)));
path5      =  1*(a*t.^2-y4);

%----------------------------------------------------%
% Compute the Differential Equations in Phase IPHASE %
%----------------------------------------------------%
if ~(iphase==3),
  path = [path1,path2,path3,path4,path5];
else
  path = [path1,path2,path3,path4];
end;  
phaseout(iphase).dynamics  = [y1dot,y2dot,y3dot,y4dot,y5dot,y6dot];
phaseout(iphase).path = path;

end

%------------------------------------%
% END: kineticBatchReactorContinuous.m
%------------------------------------%
