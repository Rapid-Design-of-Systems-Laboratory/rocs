function [p,f,g,h,k,L] = lowThrustCoe2Mee(a,ecc,inc,Ome,ome,nu)

p = a*(1-ecc^2);
f = ecc*cos(ome+Ome);
g = ecc*sin(ome+Ome);
h = tan(inc/2)*cos(Ome);
k = tan(inc/2)*sin(Ome);
L = Ome+ome+nu;

end
