%Brad Steinfeldt
%D'Souza Lunar Lander Example

function LunarLanding

close all
clc
clear all

%Initial Conditions
r0 = 1000*[500; 100; 50];        %ft
v0 =  1000*[-3; 0; 0];                %ft/s

%Target Conditions
rT = [0; 0; 0];                     %ft
vT = [0; 0; 0];                     %ft/s

%Problem Constants
g = 5.32;                           %ft/s^2
W = 0;                              %Weighting on the time-to-go

%Form IC Vector
X0 = [r0;v0];

%Setup the Integration Routine
tspan = [0 10000];
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'events',@HitGround);
[T,X,TE,XE,IE] = ode45(@LunarEOMs,tspan,X0,options,rT,vT,g,W);

%Post Process Accelerations
for i = 1 : size(T,1)
    
    x = X(i,1);
    y = X(i,2);
    z = X(i,3);
    u = X(i,4);
    v = X(i,5);
    w = X(i,6);

    rV = [x;y;z];
    vV = [u;v;w];

    a_commanded =  DSouzaAcceleration(rV,vV,rT,vT,g,W);

    ax(i) = a_commanded(1);
    ay(i) = a_commanded(2);
    az(i) = a_commanded(3);
    
    V(i) = norm(vV);
    
end
    
%Plot CR vs DR
subplot(2,2,1)
plot(X(:,1)/1000,X(:,2)/1000);
xlabel('Downrange [kft]');
ylabel('Crossrange [kft]');
grid on

%Plot h vs DR
subplot(2,2,2)
plot(X(:,1)/1000,X(:,3)/1000);
xlabel('Downrange [kft]');
ylabel('Altitude [kft]');
grid on

%Plot h vs V
subplot(2,2,[3 4])
plot(V,X(:,3)/1000);
xlabel('Velocity [ft/s]');
ylabel('Altitude [kft]');
grid on

figure

%Plot Accelerations
subplot(3,1,1)
plot(T(:,1),ax);
xlabel('Time [s]');
ylabel('x-acceleration [ft/s^2]');
grid on
subplot(3,1,2)
plot(T(:,1),ax);
xlabel('Time [s]');
ylabel('y-acceleration [ft/s^2]');
grid on
subplot(3,1,3)
plot(T(:,1),az);
xlabel('Time [s]');
ylabel('z-acceleration [ft/s^2]');
grid on

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xdot] = LunarEOMs(t,X,rT,vT,g,W)

x = X(1);
y = X(2);
z = X(3);
u = X(4);
v = X(5);
w = X(6);

rV = [x;y;z];
vV = [u;v;w];

a_commanded =  DSouzaAcceleration(rV,vV,rT,vT,g,W);

rdot = [u;v;w];
vdot = a_commanded + [0;0;g];

Xdot = [rdot;vdot];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction] = HitGround(t,X,rT,vT,g,W)

x = X(1);
y = X(2);
z = X(3);
u = X(4);
v = X(5);
w = X(6);

value = z;
direction = 0;
isterminal = 1;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a_commanded] = DSouzaAcceleration(rV,vV,rT,vT,g,W)

%Compute target relative distance and velocity
r = rV - rT;
v = vV - vT;

%Compute quartic coefficients for the time-to-go tgo^4 + a*tgo^2 + b*tgo +
%c =
a = - (2*dot(v,v))/(W+g^2/2);
b = - (12*dot(v,r))/(W+g^2/2);
c = - (18*dot(r,r))/(W+g^2/2);

TGO = roots([1 0 a b c]);

tgo = 1e30;
for i = 1 : length(TGO)
    if isreal(TGO(i)) && real(TGO(i)) > 0 && real(TGO(i)) < tgo
        tgo = TGO(i);
    end
end

a_commanded = -4*v/tgo - 6*r/tgo^2 - [0;0;g];

return