function Y=blackbox32(P,X)
%BLACKBOX32
%  Interplate very accurate solutions from MIRK12 for CW32

%   Nick Hale
%       Imperial College London 
%       June 2006

coeffs;     %compute constant coefficients

if P*(33-P)<0
    error('Problem outside permited range: 0-32\n');
end 
[XL,XR]=xlist(P);
E=epslist(P);
N=length(X);

% exist analytic solution?
switch P 
    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20}
        Y=soln(X,E,P);  
        return; 
end;

p=0;
info=csvread('CWBVP_data.txt',0,0,[0,0,0,2]);
while info(1)~=P
    p=p+info(2)+1;
    info=csvread('CWBVP_data.txt',p,0,[p,0,p,2]);
end
y=csvread('CWBVP_data.txt',p+1,0,[p+1,0,p+info(2),info(3)-1]); %mirk 12 solns
x=XL:(XR-XL)/(info(3)-1):XR;

if all(XL<=X) + all(X<=XR)<2    %check that XL<X<XR
    error('X outside interval range [%f,%f]\n',XL,XR);
end

strt=1;
if X(1)==XL                     %check if X[1]==XL
    I(1)=0; w(1)=0; strt=2;
end
for i=strt:N
    I(i)=0;
    if X(i)==XL                  %error if X[i]==XL for i!=1
        error('X[i]=XL for i!=1. Please order the X[i]');
    end
    while X(i)-x(I(i)+1)>0      %find interval containing each X[i]
        I(i)=I(i)+1;
    end
    w(i)=(X(i)-x(I(i)))./(x(I(i)+1)-x(I(i)));
end

J=1;K=1;
if I(J)==0                      %if X[J]=XL
    Y(:,J)=y(:,1);
    J=2;K=2;
end
while J<=N
    while I(K)-I(J)==0          %find X's in same interval
        K=K+1;
        if K==N+1 break; end
    end
 
    if w(K-1)==1                %if X[J]=x[j] then use exact y(:,j)
        Y(:,(K-1))=y(:,I(J)+1);
        if J<=K-2 
            Y(:,J:(K-2))=interp(x,y,E,P,w(J:(K-2))',I(J));
        end
    else
        Y(:,J:(K-1))=interp(x,y,E,P,w(J:(K-1))',I(J));
    end
    J=K;
end

% %---------------------------------------------------------- %delete this
X;
Y;


%-----------------------------------------------------------
  
function Y=interp(x,y,E,P,w,int)
XL=x(int);
XR=x(int+1);
YL=y(:,int);
YR=y(:,int+1);
FL=f(XL,YL,E,P);
FR=f(XR,YR,E,P);

%find interpolation points
[ypa7p,yma7p,ypa8p,yma8p,ypa9p,yma9p]=interppts(XL,XR,YL,YR,FL,FR,E,P);
%interpolate at w        
for i=1:length(w)
    wi=w(i);
    Awi=A(wi);
    Y(:,i)=Awi*YR+(1-Awi)*YL+     ...
       (XR-XL)*( B(wi)*FR-B(1-wi).*FL          ...
                +C(wi)*ypa7p-C(1-wi)*yma7p    ...
                +D(wi)*ypa8p-D(1-wi)*yma8p    ...
                +EE(wi)*ypa9p-EE(1-wi)*yma9p);
end



 
function [ypa7p,yma7p,ypa8p,yma8p,ypa9p,yma9p]=interppts(XL,XR,YL,YR,FL,FR,E,P)
%---------------------------------COEFFS-----------------------------------
global r3 r5 r11 r15 r33 r1311 r1653 r2185 r4054 r22857 r24795 r5294425981
global a1 a2 a3 a4 a5 a6 a7 a8 a9
global a1m a1p a2m a2p
global b1p b1m b2p b2m b3
global c1p c1m c2p c2m c3p c3m
global d1p d1m d2p d2m d3p d3m
global e1 e2 e3
global f1p f1m f2p f2m f3 f4p f4m f5
global g1p g1m g2p g2m g3 g4p g4m g5
global rr1 rr2 rr3 rr4
global s1p s1m s2p s2m s3p s3m s4p s4m
global t1p t1m t2p t2m t3p t3m t4p t4m
global u1p u1m u2p u2m u3 u4 u5p u5m u6m u6p
%--------------------------------------------------------------------------
h=XR-XL;
XLh=XL+0.5*h;

ypa1p=f(XLh+h*a1,a1p*YL+a1m*YR-h*(a2p*FL+a2m*FR),E,P);
yma1p=f(XLh-h*a1,a1m*YL+a1p*YR+h*(a2m*FL+a2p*FR),E,P);

ypa2p=f(XLh+h*a2,b1p*YL+b1m*YR-h*(b2p*FL+b2m*FR+b3*(yma1p-ypa1p)),E,P);
yma2p=f(XLh-h*a2,b1m*YL+b1p*YR+h*(b2m*FL+b2p*FR-b3*(yma1p-ypa1p)),E,P);

ypa3p=f(XLh+h*a3,c1p*YL+c1m*YR-h*(c2p*FL+c2m*FR+c3p*yma2p+c3m*ypa2p),E,P);
yma3p=f(XLh-h*a3,c1m*YL+c1p*YR+h*(c2m*FL+c2p*FR+c3m*yma2p+c3p*ypa2p),E,P);

ypa4p=f(XLh+h*a4,d1p*YL+d1m*YR-h*(d2p*FL+d2m*FR+d3p*yma3p+d3m*ypa3p),E,P);
yma4p=f(XLh-h*a4,d1m*YL+d1p*YR+h*(d2m*FL+d2p*FR+d3m*yma3p+d3p*ypa3p),E,P);

yp05p=f(XLh,0.5*(YL+YR)+h*(e1*(FL-FR)+e2*(yma3p-ypa3p)+e3*(yma4p-ypa4p)),E,P);

ypa5p=f(XLh+h*a5,f1p*YL+f1m*YR-h*(f2p*FL+f2m*FR+f3*(ypa3p-yma3p)+f4p*yma4p+f4m*ypa4p+f5*yp05p),E,P);
yma5p=f(XLh-h*a5,f1m*YL+f1p*YR+h*(f2m*FL+f2p*FR-f3*(ypa3p-yma3p)+f4m*yma4p+f4p*ypa4p+f5*yp05p),E,P);

ypa6p=f(XLh+h*a6,g1p*YL+g1m*YR-h*(g2p*FL+g2m*FR+g3*(ypa3p-yma3p)+g4p*yma4p+g4m*ypa4p+g5*yp05p),E,P);
yma6p=f(XLh-h*a6,g1m*YL+g1p*YR+h*(g2m*FL+g2p*FR-g3*(ypa3p-yma3p)+g4m*yma4p+g4p*ypa4p+g5*yp05p),E,P);

yp05p=f(XLh,0.5*(YL+YR)+h*(rr1*(FL-FR)+rr2*(yma4p-ypa4p)+rr3*(yma5p-ypa5p)+rr4*(yma6p-ypa6p)),E,P);

ypa7p=f(XLh+h*a7,s1p*YL+s1m*YR-h*(s2p*FL+s2m*FR+s3p*yma5p+s3m*ypa5p+s4p*yma6p+s4m*ypa6p),E,P);
yma7p=f(XLh-h*a7,s1m*YL+s1p*YR+h*(s2m*FL+s2p*FR+s3m*yma5p+s3p*ypa5p+s4m*yma6p+s4p*ypa6p),E,P);

ypa8p=f(XLh+h*a8,t1p*YL+t1m*YR-h*(t2p*FL+t2m*FR+t3p*yma5p+t3m*ypa5p+t4p*yma6p+t4m*ypa6p),E,P);
yma8p=f(XLh-h*a8,t1m*YL+t1p*YR+h*(t2m*FL+t2p*FR+t3m*yma5p+t3p*ypa5p+t4m*yma6p+t4p*ypa6p),E,P);

ypa9p=f(XLh+h*a9,u1p*YL+u1m*YR-h*(u2p*FL+u2m*FR+u3*(ypa5p-yma5p)+u4*(ypa6p-yma6p)+ ...
    u5p*yma7p+u5m*ypa7p+u6p*yma8p+u6m*ypa8p),E,P);
yma9p=f(XLh-h*a9,u1m*YL+u1p*YR+h*(u2m*FL+u2p*FR-u3*(ypa5p-yma5p)-u4*(ypa6p-yma6p)+ ...
    u5m*yma7p+u5p*ypa7p+u6m*yma8p+u6p*ypa8p),E,P);

%----------------------------Interpolation functions----------------------------------------    

function coeffs = A(w)
    coeffs = -w.^2.*(-210+2590*w-14070*w.^2+41664*w.^3-71610*w.^4+71280*w.^5-38115*w.^6+8470*w.^7);
function coeffs = B(w)
    coeffs = 1/24*w.^2.*(w-1).*(8833*w.^6-30371*w.^5+42097*w.^4-30008*w.^3+11620*w.^2-2354*w+207);
function coeffs = C(w)
    global r3 r11 r33 a8
    coeffs = 1/95568*( 1445466*w.^5+(91476*a8*r11-176418*a8*r33-3613665)*w.^4  ...
                   +(-182952*a8*r11+352836*a8*r33+8712*r3+3441240)*w.^3    ...
                   +(-13068*r3-265650*a8*r33-1548195+157212*a8*r11)*w.^2   ...
                   +(89232*a8*r33+334158-65736*a8*r11+8316*r3)*w            ...
                   -29502-1980*r3+11880*a8*r11-12672*a8*r33).*(w-1).^2.*w.^2*(148+r3);
function coeffs = D(w)
    global r3 a8
    coeffs = 1/95568*w.^2.*(w-1).^2.*(r3-148).* ...
                (-1445466*w.^5 ...
                 +(3613665+1221858*a8-431244*a8*r3)*w.^4 ...
                 +(-2443716*a8+862488*a8*r3+8712*r3-3441240)*w.^3 ...
                 +(-13068*r3-590964*a8*r3+1548195+1762002*a8)*w.^2 ...
                 +(159720*a8*r3-334158+8316*r3-540144*a8)*w ...
                 +29502-15048*a8*r3-1980*r3+66528*a8);
function coeffs = EE(w)
    global r33
    coeffs = 121/48*w.^2.*(w-1).^2.*(-242*w.^5+(33*r33+605)*w.^4-(66*r33+528)*w.^3 ...
                                   +(49*r33+187)*w.^2-(16*r33+22)*w+2*r33);
                               
%-----------------------------------------------------------------------------------------

function dydx = f(x,y,E,n)
[ode,ignored,ignored2] = feval(@problemlist,n);
dydx = feval(ode,x,y,E);

function ssoln = soln(x,E,n)
[ignored,ignored2,sol] = feval(@problemlist,n);
ssoln = feval(sol,x,E);

%-----------------------calculate constant coeffs-----------------------------------------
function coeffs
global r3 r5 r11 r15 r33 r1311 r1653 r2185 r4054 r22857 r24795 r5294425981
r3=sqrt(3.0);
r5=sqrt(5.0);
r11=sqrt(11.0);
r15=r3*r5;
r33=r11*r3;
r1311=sqrt(1311.0);
r1653=sqrt(1653.0);
r2185=sqrt(2185.0);
r4054=sqrt(4054.0);
r22857=sqrt(22857.0);
r24795=sqrt(24795.0);
r5294425981=sqrt(5294425981.0);

global a1 a2 a3 a4 a5 a6 a7 a8 a9
a1=1/4.0;
a2=r5294425981/212914.0;
a3=r2185/114.0;
a4=r1653/114.0;
a5=sqrt(5445/(15+2*r15))/66.0;
a6=sqrt(495+66*r15)/66.0;
a7=sqrt(297-132*r3)/66.0;
a8=sqrt(297+132*r3)/66.0;
a9=r33/22.0;

global a1m a1p a2m a2p
a1m = 27/32.0;
a1p = 5/32.0;
a2m = 9/64.0;
a2p = -3/64.0;

global b1p b1m b2p b2m b3
b1p = 1/2.0-134819/22666185698.0*r5294425981;
b1m = 1/2.0+134819/22666185698.0*r5294425981;
b2p = 14181/22666185698.0*r5294425981-973398021/22666185698.0;
b2m = 14181/22666185698.0*r5294425981+973398021/22666185698.0;
b3  = -536268696/11333092849.0;

global c1p c1m c2p c2m c3p c3m
c1p = 1/2.0-154085/13678632.0*r2185;
c1m = 1/2.0+154085/13678632.0*r2185;
c2p = 71731079/122511587616.0*r2185-618856/21824559.0;
c2m = 71731079/122511587616.0*r2185+618856/21824559.0;
c3p = -274547/1085400792747.0*r5294425981+1538286841/2327720164704.0*r2185;
c3m = 274547/1085400792747.0*r5294425981+1538286841/2327720164704.0*r2185;

global d1p d1m d2p d2m d3p d3m
d1p = 1/2.0-8053/656298.0*r1653;
d1m = 1/2.0+8053/656298.0*r1653;
d2p = 1001/2625192.0*r1653-7/456.0;
d2m = 1001/2625192.0*r1653+7/456.0;
d3p = -21/17480.0*r2185+21/15352.0*r1653;
d3m = +21/17480.0*r2185+21/15352.0*r1653;

global e1 e2 e3
e1 = 827/12544.0;
e2 = -1539/288512.0*r2185;
e3 = 57/6272.0*r1653;

global f1p f1m f2p f2m f3 f4p f4m f5
f1p = (15972-5*a5*(13302+301*r15))/31944.0;
f1m = (15972+5*a5*(13302+301*r15))/31944.0;
f2p = -(3191/149072+2225/1565256*r15)+(40085/447216+1979/894432*r15)*a5;
f2m = 3191/149072+2225/1565256*r15+(40085/447216+1979/894432*r15)*a5;
f3  = 171/24000592*r1311*(7*r15-124);
f4p = -(437/308792*r1653+6745/15130808*r24795)+(1543275/4323088+267501/8646176*r15)*a5;
f4m = 437/308792*r1653+6745/15130808*r24795+(1543275/4323088+267501/8646176*r15)*a5;
f5  = 2*a5*(994-101*r15)/10527;

global g1p g1m g2p g2m g3 g4p g4m g5
g1p = 1/2-11085*a6/5324+1505/31944*r15*a6;
g1m = (15972+5*a6*(13302-301*r15))/31944.0;
g2p = -3191/149072+2225/1565256*r15+(40085/447216-1979/894432*r15)*a6;
g2m = 3191/149072-2225/1565256*r15+(40085/447216-1979/894432*r15)*a6;
g3  = 171/24000592*r1311*(7*r15+124);
g4p = -437/308792*r1653+6745/15130808*r24795+(1543275/4323088-267501/8646176*r15)*a6;
g4m = 437/308792*r1653-6745/15130808*r24795+(1543275/4323088-267501/8646176*r15)*a6;
g5  = 2*a6*(994+101*r15)/10527;

global rr1 rr2 rr3 rr4
rr1 = -121/7168;
rr2 = -1083/207872*r1653;
rr3 = 3*a5*(7*r15+124)/512;
rr4 = -3*a6*(7*r15-124)/512;

global s1p s1m s2p s2m s3p s3m s4p s4m
s1p = 1/43454004.0*(3627+874*r3)*(7254-1748*r3-29927*a7);
s1m = 1/43454004.0*(3627+874*r3)*(7254-1748*r3+29927*a7);
s2p = 1/618552.0*(171+34*r3)*(213*a7+12*r3-52);
s2m = 1/618552.0*(171+34*r3)*(213*a7-12*r3+52);
s3p = 1/19021200.0*(3340+9*r15+362*r5+1400*r3)*(1965*a7-56*r15*a5-1360*a5 ...
        -32*r5*a5+240*r3*a5);
s3m = 1/19021200.0*(3340+9*r15+362*r5+1400*r3)*(1965*a7+56*r15*a5+1360*a5 ...
        +32*r5*a5-240*r3*a5); 
s4p = 1/19021200.0*(3340+1400*r3-362*r5-9*r15)*(1965*a7-1360*a6+240*r3*a6 ...
        +32*r5*a6+56*r15*a6); 
s4m = 1/19021200.0*(3340+1400*r3-362*r5-9*r15)*(1965*a7+1360*a6-240*r3*a6 ...
        -32*r5*a6-56*r15*a6);         

global t1p t1m t2p t2m t3p t3m t4p t4m
t1p = 1/43454004.0*(874*r3-3627)*(-7254+29927*a8-1748*r3);
t1m = -1/43454004.0*(874*r3-3627)*(29927*a8+7254+1748*r3);
t2p = -1/618552.0*(-171+34*r3)*(-52-12*r3+213*a8);
t2m = -1/618552.0*(-171+34*r3)*(52+12*r3+213*a8);
t3p = 1/19021200.0*(-1400*r3+9*r5*r3+3340-362*r5)*(1965*a8+32*r5*a5 ...
        -240*r3*a5-1360*a5-56*a5*r15);
t3m = 1/19021200.0*(-1400*r3+9*r15+3340-362*r5)*(1965*a8+1360*a5 ...
        -32*r5*a5+240*r3*a5+56*a5*r15);
t4p = -1/19021200.0*(-3340+9*r5*r3+1400*r3-362*r5)*(-1360*a6 ...
        +1965*a8-240*r3*a6+56*r15*a6-32*r5*a6);
t4m = -1/19021200.0*(-3340+9*r5*r3+1400*r3-362*r5)*(240*r3*a6 ...
        +32*r5*a6+1360*a6+1965*a8-56*r15*a6);

global u1p u1m u2p u2m u3 u4 u5p u5m u6m u6p
u1m = -51/2662.0*r11*r3+1/2.0;
u1p = 51/2662.0*r11*r3+1/2.0;
u2m = 83/6655.0-1/1331.0*r11*r3;
u2p = -83/6655.0-1/1331.0*r11*r3;
u3 = 4/33275.0*(1416+37*r15)*a5;
u4 = 4/33275.0*r11*(122*r15+207)*a5;
u5m = 1/945010.0*r11*(-171+34*r3)*(120*r3+165+142*a8);
u5p = 1/945010.0*r11*(-171+34*r3)*(120*r3+165-142*a8);
u6m = 1/10395110.0*(-378+377*r3)*(15*r11*r3-420*r11+1562*a8);
u6p = 1/10395110.0*(-378+377*r3)*(-1562*a8+15*r33-420*r11);
