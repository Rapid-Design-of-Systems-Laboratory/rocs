clear all; close all; clc;
%Initialization of PSO parameters
wmax=0.9;
wmin=0.4;
itmax=200; %Maximum iteration number
c1=1.4;
c2=1.4;

for iter=1:itmax
W(iter)=wmax-((wmax-wmin)/itmax)*iter;
end

%**********************************************************

%Initialization of positions of agents
% agents are initialized between -5,+5 randomly
a=0;  
b=5;
N=20;
D=2;

x=a+(b-a)*rand(N,D,1);

%Initialization of velocities of agents
%Between -5 , +5, (which can also be started from zero)

m=0;
n=1;
V=m+(n-m)*rand(N,D,1);

%**********************************************************
%Function to be minimized. Scheaffer F6
%x(i,2,1)=x2
%x(i,1,1)=x1
%f=100*(x2-x1^2)^2+(x1-1)^2
for i=1:N;
F(i,1,1)=100*(x(i,2,1)-x(i,1,1)^2)^2+(x(i,1,1)-1)^2;
end
%**********************************************************
[C,I]=min(abs(F(:,1,1)));

B(1,1,1)=C;
XX(1,1,1)=I;
gbest(1,1,1)=x(I,1,1);
gbest(1,2,1)=x(I,2,1);

%********************************************************
%Matrix composed of gbest vector 
for p=1:N
    for r=1:D  
    G(p,r,1)=gbest(1,r,1);
end
  end
 
  Fbest(1,1,1)=100*(G(1,2,1)-G(1,1,1)^2)^2+(G(1,1,1)-1)^2;

  for i=1:N;
       pbest(i,:,1)=x(i,:,1);
    end
    

V(:,:,2)=W(1)*V(:,:,1)+c1*rand*(pbest(:,:,1)-x(:,:,1))+c2*rand*(G(:,:,1)-x(:,:,1));

x(:,:,2)=x(:,:,1)+V(:,:,2);

Fb(1,1,1)=100*(gbest(1,2,1)-gbest(1,1,1)^2)^2+(gbest(1,1,1)-1)^2;

%******************************************************
for j=2:itmax-1
   
% Calculation of new positions

    
for i=1:N;
    F(i,1,j)=100*(x(i,2,j)-x(i,1,j)^2)^2+(x(i,1,j)-1)^2;
end


[C,I]=min(abs(F(:,:,j)));

B(1,1,j)=C;

gbest(1,1,j)=x(I,1,j);
gbest(1,2,j)=x(I,2,j);
        
Fb(1,1,j)=100*(gbest(1,2,j)-gbest(1,1,j)^2)^2+(gbest(1,1,j)-1)^2;

[C,I]=min(Fb(1,1,:));

if Fb(1,1,j)<=C
    gbest(1,1,j)=gbest(1,1,j);
    gbest(1,2,j)=gbest(1,2,j);
else
    gbest(1,1,j)=gbest(1,1,I);
    gbest(1,2,j)=gbest(1,2,I);
end
   
    
%Matrix composed of gbest vector 
for p=1:N
    for r=1:D  
    G(p,r,j)=gbest(1,r,j);
end
  end

 Fbest(1,1,j)=100*(G(1,2,j)-G(1,1,j)^2)^2+(G(1,1,j)-1)^2;

 
  for i=1:N;
      [C,I]=min(F(i,1,:));
    if F(i,1,j)<=C
       pbest(i,:,j)=x(i,:,j);
   else
        pbest(i,:,j)=x(i,:,I);
    end
end

V(:,:,j+1)=W(j)*V(:,:,j)+c1*rand*(pbest(:,:,j)-x(:,:,j))+c2*rand*(G(:,:,j)-x(:,:,j));

x(:,:,j+1)=x(:,:,j)+V(:,:,j+1);

end


plot(x(:,1,200),x(:,2,200),'o')

