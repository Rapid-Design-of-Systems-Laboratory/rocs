function f=functn(nopt,x)
%
%  This is the Shekel Function
%  Bound: X(i)=[0,10], j=1,2,3,4
%  Global Optimum:-10.5364098252,(4,4,4,4)
%
%  Data for Skekel function coefficients (n=4, m=10)
      a1=[4.,1.,8.,6.,3.,2.,5.,8.,6.,7.; ...
          4.,1.,8.,6.,7.,9.,5.,1.,2.,3.6; ...
          4.,1.,8.,6.,3.,2.,3.,8.,6.,7.; ...
          4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6];
      c1 =[.1,.2,.2,.4,.4,.6,.3,.7,.5,.5];
%
f = 0.0;
for i = 1:10
    u = 0.0;
    for j = 1:nopt
        u = u + (x(j) - a1(j,i))^2;
    end;
    u = 1.0 / (u + c1(i));
    f = f - u;
end;

return
