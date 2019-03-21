function [a,e,i,Ome,ome,nu] = lowThrustMee2Coe(p,f,g,h,k,L)

a = p./(1-f.*f-g.*g);
e = sqrt(f.*f+g.*g);
i = atan2(2.*sqrt(h.*h+k.*k),1-h.*h-k.*k);
ome = atan2(g.*h-f.*k,f.*h+g.*k);
Ome = atan2(k,h);
nu = L-Ome-ome;

end
