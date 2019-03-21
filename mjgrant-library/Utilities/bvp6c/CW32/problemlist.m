function [dydx,res,ans] = problemlist(n)
switch n
    case 1
        dydx = @ode1;
        res = @bc1;
        ans = @soln1;
    case 2 
        dydx = @ode2;
        res = @bc2;
        ans = @soln2;
    case 3
        dydx = @ode3;
        res = @bc3;
        ans = @soln3;
    case 4
        dydx = @ode4;
        res = @bc4;
        ans = @soln4;
    case 5
        dydx = @ode5;        
        res = @bc5;
        ans = @soln5;
    case 6
        dydx = @ode6;
        res = @bc6;
        ans = @soln6;
    case 7
        dydx = @ode7;
        res = @bc7;
        ans = @soln7;
    case 8
        dydx = @ode8;        
        res = @bc8;
        ans = @soln8;
    case 9
        dydx = @ode9;
        res = @bc9;
        ans = @soln9;
    case 10
        dydx = @ode10;
        res = @bc10;
        ans = @soln10;
    case 11
        dydx = @ode11;        
        res = @bc11;
        ans = @soln11;
    case 12
        dydx = @ode12;
        res = @bc12;
        ans = @soln12;
    case 13
        dydx = @ode13;
        res = @bc13;
        ans = @soln13;
    case 14
        dydx = @ode14;        
        res = @bc14;
        ans = @soln14;
    case 15
        dydx = @ode15;
        res = @bc15;
        ans = @soln15;
    case 16
        dydx = @ode16;  
        res = @bc16;
        ans = @soln16;
    case 17
        dydx = @ode17;  
        res = @bc17;
        ans = @soln17;
    case 18
        dydx = @ode18;  
        res = @bc18;
        ans = @soln18;
    case 19
        dydx = @ode19;  
        res = @bc19;
        ans = @soln19;
    case 20
        dydx = @ode20;  
        res = @bc20;
        ans = @soln20;
    case 21
        dydx = @ode21;  
        res = @bc21;
        ans = @soln21;
    case 22
        dydx = @ode22;  
        res = @bc22;
        ans = @soln22;
    case 23
        dydx = @ode23;  
        res = @bc23;
        ans = @soln23;
    case 24
        dydx = @ode24;  
        res = @bc24;
        ans = @soln24;
    case 25
        dydx = @ode25;
        res = @bc25;
        ans = @soln25;
    case 26
        dydx = @ode26;          
        res = @bc26;
        ans = @soln26;
    case 27
        dydx = @ode27;  
        res = @bc27;
        ans = @soln27;
    case 28
        dydx = @ode28;  
        res = @bc28;
        ans = @soln28;
    case 29
        dydx = @ode29;  
        res = @bc29;
        ans = @soln29;
    case 30
        dydx = @ode30;  
        res = @bc30;
        ans = @soln30;
    case 31
        dydx = @ode31;  
        res = @bc31;
        ans = @soln31;
    case 32
        dydx = @ode32;          
        res = @bc32;
        ans = @soln32;
end

%----------------------------------------   Problem 1 
function dydx = ode1(x,y,E)

  dydx = [ y(2) 
           y(1)/E];

function res = bc1(ya,yb,E)
  res = [ ya(1) - 1
          yb(1)];
      
function ans = soln1(x,E)
  sE = sqrt(E);
  ans(1,:) = (exp(-x/sE)-exp((x-2)/sE)) /(1-exp(-2/sE));
  ans(2,:) = (exp(-x/sE)+exp((x-2)/sE))/(sE*(exp(-2/sE)-1));
  
%----------------------------------------   Problem 2 
  
function dydx = ode2(x,y,E)
  dydx = [ y(2) 
           y(2)/E];

function res = bc2(ya,yb,E)
  res = [ ya(1) - 1
          yb(1)];
      
function ans = soln2(x,E)
  ans(1,:) = (1-exp((x-1)/E))/(1-exp(-1/E));   
  ans(2,:) = -exp((x-1)/E)/(E*(1-exp(-1/E)));
  
%-----------------------------------------   Problem 3

function dydx = ode3(x,y,E)
  R = -(1.0+E*pi^2)*cos(pi*x)-(2.0+cos(pi*x))*pi*sin(pi*x);
  dydx = [ y(2) 
           (R+y(1)-(2.0+cos(pi*x))*y(2))/E];

function res = bc3(ya,yb,E)
  res = [ ya(1) + 1
          yb(1) + 1];
      
function ans = soln3(x,E)
  ans(1,:) = cos(pi*x);   
  ans(2,:) = -pi*sin(pi*x);
  
%----------------------------------------   Problem  4

function dydx = ode4(x,y,E)
  dydx = [ y(2) 
           ((1+E)*y(1)-y(2))/E];

function res = bc4(ya,yb,E)
  res = [ ya(1) - (1+exp(-2))
          yb(1) - (1+exp(-2*(1+E)/E))];
      
function ans = soln4(x,E)
  ans(1,:) = exp(x-1)+exp(-(1+E)*(1+x)/E);   
  ans(2,:) = exp(x-1.0)-((1+E)*exp(-(1+E)*(1+x)/E))/E;
  
%----------------------------------------   Problem  5
function dydx = ode5(x,y,E)
  dydx = [ y(2) 
           (y(1)+x*y(2)-(1+E*pi^2)*cos(pi*x)+pi*x*sin(pi*x))/E];

function res = bc5(ya,yb,E)
  res = [ ya(1) + 1
          yb(1) + 1];
      
function ans = soln5(x,E)
  ans(1,:) = cos(pi*x);   
  ans(2,:) = -pi*sin(pi*x);
  
%----------------------------------------   Problem  6

function dydx = ode6(x,y,E)
  dydx = [ y(2) 
           (-x*y(2)-E*pi^2*cos(pi*x)-pi*x*sin(pi*x))/E];

function res = bc6(ya,yb,E)
  res = [ ya(1) + 2
          yb(1)];
      
function ans = soln6(x,E)
  ans(1,:) = cos(pi*x)+erf(x/sqrt(2*E))/erf(1/sqrt(2*E));   
  ans(2,:) = -sin(pi*x)*pi+sqrt(2/(pi*E))*exp(-x.^2/(2*E))/erf(1/sqrt(2*E));
  
function dydx = ode(x,y,E)
  dydx = [ y(2) 
           (y(1)-x*y(2)-(1+E*pi^2)*cos(pi*x)-pi*x*sin(pi*x))/E];

%----------------------------------------   Problem   7

function dydx = ode7(x,y,E)
  dydx = [ y(2) 
           (y(1)-x*y(2)-(1+E*pi^2)*cos(pi*x)-pi*x*sin(pi*x))/E];

function res = bc7(ya,yb,E)
  res = [ ya(1) + 1
          yb(1) - 1];
      
function ans = soln7(x,E)
  ans(1,:) = cos(pi*x)+x+ (x.*erf(x/sqrt(2*E))+sqrt(2*E/pi)*exp(-x.^2/(2*E))) ...
      / (erf(1/sqrt(2*E)) + sqrt(2*E/pi)*exp(-1/(2*E)));  
  ans(2,:) = 1-sin(pi*x)*pi+erf(x/sqrt(2*E)) ...
      / (erf(1/sqrt(2*E)) + sqrt(2*E/pi)*exp(-1/(2*E)));
  
%----------------------------------------   Problem  8

function dydx = ode8(x,y,E)
  dydx = [ y(2) 
           -y(2)/E];

function res = bc8(ya,yb,E)
  res = [ ya(1) - 1
          yb(1) - 2];
      
function ans = soln8(x,E)
  ans(1,:) = (2-exp(-1/E)-exp(-x/E))/(1-exp(-1/E)); 
  ans(2,:) = exp(-x/E)/(1-exp(-1/E))/E;
  
%----------------------------------------   Problem  9

function dydx = ode9(x,y,E)
  dydx = [ y(2) 
           -(2*y(1)+4*x*y(2))./(E+x.^2)];

function res = bc9(ya,yb,E)
  res = [ ya(1) - 1/(1+E)
          yb(1) - 1/(1+E)];
      
function ans = soln9(x,E)
  ans(1,:) = 1./(E+x.^2);
  ans(2,:) = -2*x./(E+x.^2).^2;
  
%----------------------------------------   Problem  10

function dydx = ode10(x,y,E)
  dydx = [ y(2) 
           -x*y(2)/E];

function res = bc10(ya,yb,E)
  res = [ ya(1)
          yb(1) - 2];
      
function ans = soln10(x,E)
  ans(1,:) = 1+erf(x/sqrt(2*E))/erf(1/sqrt(2*E));
  ans(2,:) = sqrt(2/(E*pi))*exp(-x.^2/(2*E))/erf(1/sqrt(2*E));
  
%----------------------------------------   Problem 11   

function dydx = ode11(x,y,E)
  dydx = [ y(2) 
           (y(1)-(E*pi^2+1)*cos(pi*x))/E];

function res = bc11(ya,yb,E)
  res = [ ya(1) + 1
          yb(1) + 1];
      
function ans = soln11(x,E)
  ans(1,:) = cos(pi*x); 
  ans(2,:) = -pi*sin(pi*x);
  
%----------------------------------------   Problem   12

function dydx = ode12(x,y,E)
  dydx = [ y(2) 
           (y(1)-(E*pi^2+1)*cos(pi*x))/E];

function res = bc12(ya,yb,E)
  res = [ ya(1) + 1
          yb(1)];
      
function ans = soln12(x,E)
    sE = sqrt(E);
    ans(1,:) = cos(pi*x)+(exp((x+1)/sE)-exp((-x-1)/sE))/ ...
               (exp(2/sE)-exp(-2/sE));
    ans(2,:) = -pi*sin(pi*x)+(exp((x+1)/sE)+exp((-x-1)/sE))/ ...
               (exp(2/sE)-exp(-2/sE))/sE;
           
         
%----------------------------------------   Problem  13

function dydx = ode13(x,y,E)
  dydx = [ y(2) 
           (y(1)-(E*pi^2+1)*cos(pi*x))/E];

function res = bc13(ya,yb,E)
  res = [ ya(1)
          yb(1) + 1];
      
function ans = soln13(x,E)
  
  sE = sqrt(E);
  A = exp(-1/sE)/(exp(2/sE)-exp(-2/sE));
  B = exp(1/sE)/(exp(2/sE)-exp(-2/sE));
        
  ans(1,:) = -A*exp(x/sE)+B*exp(-x/sE)+cos(pi*x);
  ans(2,:) = (-A*exp(x/sE)-B*exp(-x/sE))/sE-pi*sin(pi*x);
  
%----------------------------------------   Problem  14
function dydx = ode14(x,y,E)
  dydx = [ y(2) 
           (y(1)-(E*pi^2+1)*cos(pi*x))/E];

function res = bc14(ya,yb,E)
  res = [ ya(1)
          yb(1)];
      
function ans = soln14(x,E)
  sE = sqrt(E);  
  A = exp(1/sE)+exp(-1/sE);
                                    
  ans(1,:) = (exp(x/sE)+exp(-x/sE))/A + cos(pi*x);
  ans(2,:) = (exp(x/sE)-exp(-x/sE))/(sE*A) - pi*sin(pi*x);
  
%----------------------------------------   Problem   15

function dydx = ode15(x,y,E)
  dydx = [ y(2) 
           x*y(1)/E];

function res = bc15(ya,yb,E)
  res = [ ya(1) - 1
          yb(1) - 1];
      
function ans = soln15(x,E)
  K=E^(-1/3);
  AAp=airy(K);
  AAm=airy(-K);
  ABp=airy(2,K);
  ABm=airy(2,-K);

  B=(-AAm+AAp)/(-AAm*ABp+AAp*ABm);
  A=(ABm-ABp)/(-AAm*ABp+AAp*ABm);
  ans(1,:) = real(A*airy(x*K)+B*airy(2,x*K));
  ans(2,:) = K*real(A*airy(1,x*K)+B*airy(3,x*K));

 
         
%----------------------------------------   Problem   16

function dydx = ode16(x,y,E)
  dydx = [ y(2) 
           -pi^2*y(1)/(4*E^2)];

function res = bc16(ya,yb,E)
  res = [ ya(1)
          yb(1) - sin(pi/(2*E))];
      
function ans = soln16(x,E)
  ans(1,:) = sin(pi*x/(2*E));
  ans(2,:) = cos(pi*x/(2.0*E))*pi/(2*E);
  
%----------------------------------------   Problem  17

function dydx = ode17(x,y,E)
  dydx = [ y(2) 
           -3*E*y(1)/(E+x.^2).^2];

function res = bc17(ya,yb,E)
  res = [ ya(1) + 0.1/sqrt(E+0.01)
          yb(1) - 0.1/sqrt(E+0.01)];
      
function ans = soln17(x,E)
  ans(1,:) = x./sqrt(E+x.^2);
  ans(2,:) = E./((E+x.^2).^(3.0/2.0));
  
%----------------------------------------   Problem   18

function dydx = ode18(x,y,E)
  dydx = [ y(2) 
           -y(2)/E];

function res = bc18(ya,yb,E)
  res = [ ya(1) - 1
          yb(1) - exp(-1/E)];
      
function ans = soln18(x,E)
  ans(1,:) = exp(-x/E);
  ans(2,:) = -exp(-x/E)/E;
  
%----------------------------------------   Problem   19

function dydx = ode19(x,y,E)
  dydx = [ y(2) 
           (pi/2*sin(pi*x/2)*exp(2*y(1))-exp(y(1))*y(2))/E];

function res = bc19(ya,yb,E)
  res = [ ya(1)
          yb(1)];
      
function ans = soln19(x,E)
  %approximation of true solution for initial guess   !!! 
  ans(1,:) = -log((1+cos(pi*x/2.0)).*(1-exp(-x/(2*E))/2.0));
  ans(2,:) = ((sin(pi*x/2.0)*pi*E+1+cos(pi*x/2.0)).*exp(-x/(2*E)) ...
               -2*sin(pi*x/2.0)*pi*E)./ ...
               (2*E*(1+cos(pi*x/2.0)).*(-2.0+exp(-x/(2*E)))); 
           
%----------------------------------------   Problem  20

function dydx = ode20(x,y,E)
  dydx = [ y(2) 
           (1-y(2)^2)/E];

function res = bc20(ya,yb,E)
  res = [ ya(1) - (1+E*log(cosh(-0.745/E)))
          yb(1) - (1+E*log(cosh(0.255/E)))];
      
function ans = soln20(x,E)
  ans(1,:) = 1+E*log(cosh((x-0.745)/E));
  ans(2,:) = tanh((x-0.745)/E) ;
  
%----------------------------------------   Problem  21

function dydx = ode21(x,y,E)
  dydx = [ y(2) 
           (y(1)+y(1)^2-exp(-2*x/sqrt(E)))/E];

function res = bc21(ya,yb,E)
  res = [ ya(1) - 1
          yb(1) - exp(-1/sqrt(E))];
      
function ans = soln21(x,E)
  sE = sqrt(E);  
  ans(1,:) = exp(-x/sE);
  ans(2,:) = -exp(-x/sE)/sE;
  
%----------------------------------------   Problem   22

function dydx = ode22(x,y,E)
  dydx = [ y(2) 
           -(y(2)+y(1)^2)/E];

function res = bc22(ya,yb,E)
  res = [ ya(1)
          yb(1) - 1/2];
      
function ans = soln22(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = x;
  ans(2,:) = 1;
  
  %----------------------------------------   Problem   23

function dydx = ode23(x,y,E)
  dydx = [ y(2) 
           E*sinh(E*y(1))];

function res = bc23(ya,yb,E)
  res = [ ya(1)
          yb(1)-1];
      
function ans = soln23(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = x;
  ans(2,:) = 1;
  
  %----------------------------------------   Problem   24

function dydx = ode24(x,y,E)
        lambda = 1.4;
        mP = E*(1+x.^2).*y(1);
        Q = -((1+lambda)/2.0-2.0*E*x).*y(1);
        R = y(2)/y(1)+2*x./(1+x.^2).*(1-(lambda-1)/2.0*y(1).^2);

        dydx = [ y(2) 
                (-Q*y(2)-R)/mP];

function res = bc24(ya,yb,E)
  res = [ ya(1)-0.9129
          yb(1)-0.375];
      
function ans = soln24(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = 0.9129-0.5379*x;
  ans(2,:) = -0.5379;

  %----------------------------------------   Problem  25 

function dydx = ode25(x,y,E)
  dydx = [ y(2) 
           y(1)./E*(1-y(2))];

function res = bc25(ya,yb,E)
  res = [ ya(1)+1.0/3.0
          yb(1)-1.0/3.0];
      
function ans = soln25(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = (x-0.5)*2.0/3.0;
  ans(2,:) = 2.0/3.0;
  
%----------------------------------------   Problem   26

function dydx = ode26(x,y,E)
  dydx = [ y(2) 
           y(1)./E*(1-y(2))];

function res = bc26(ya,yb,E)
  res = [ ya(1)-1.0
          yb(1)+1.0/3.0];
      
function ans = soln26(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = 1-4.0*x/3.0;
  ans(2,:) = -4.0/3.0;
  
%----------------------------------------   Problem   27

function dydx = ode27(x,y,E)
  dydx = [ y(2) 
           y(1)./E*(1-y(2))];

function res = bc27(ya,yb,E)
  res = [ ya(1)-1.0
          yb(1)-1.0/3.0];
      
function ans = soln27(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = 1-2.0*x/3.0;
  ans(2,:) = -2.0/3.0;
  
%----------------------------------------   Problem   28

function dydx = ode28(x,y,E)
  dydx = [ y(2) 
           y(1)*(1-y(2))/E];

function res = bc28(ya,yb,E)
  res = [ ya(1)-1.0
          yb(1)-3.0/2.0];
      
function ans = soln28(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = 1+x/2.0;
  ans(2,:) = 1.0/2.0;
  
%----------------------------------------   Problem   29

function dydx = ode29(x,y,E)
  dydx = [ y(2) 
           y(1)./E*(1-y(2))];

function res = bc29(ya,yb,E)
  res = [ ya(1)
          yb(1)-3.0/2.0];
      
function ans = soln29(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = 3.0*x/2.0;
  ans(2,:) = 3.0/2.0;
  
%----------------------------------------   Problem   30

function dydx = ode30(x,y,E)
  dydx = [ y(2) 
           y(1)./E*(1-y(2))];

function res = bc30(ya,yb,E)
  res = [ ya(1)+7.0/6.0
          yb(1)-3.0/2.0];
      
function ans = soln30(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = (16.0*x-7.0)/6.0;
  ans(2,:) = 8.0/3.0;
  
%----------------------------------------   Problem   31

function dydx = ode31(x,y,E)
  T = 1.0/cos(y(2))+E*y(4)*tan(y(2));
  
  dydx = [ sin(y(2))
           y(3)
           -y(4)/E
           ((y(1)-1)*cos(y(2))-T*y(3))/E];

function res = bc31(ya,yb,E)
  res = [ ya(1)
          yb(1)
          ya(3)
          yb(3)];
      
function ans = soln31(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = -x.*(x-1.0)/2.0;
  ans(2,:) = (1.0/2.0-x)*4.0/5.0;
  ans(3,:) = cos(1.5)./cos(3*x-1.5)-1.0;
  ans(4,:) = (0.5-x)*1.6;
  
%----------------------------------------   Problem   32

function dydx = ode32(x,y,E)
  dydx = [ y(3)
           y(4)
           y(2)
           E*(y(3)*y(2)-y(1)*y(4))];

function res = bc32(ya,yb,E)
  res = [ ya(1)
          ya(3)
          yb(1)-1.0
          yb(3)];
      
function ans = soln32(x,E)
%approximation of true solution for initial guess   !!! 
  ans(1,:) = x;
  ans(2,:) = 1.0;
  ans(3,:) = x.*(x-1);
  ans(4,:) = 1.0;