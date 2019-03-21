%-------------------------------------------------------------------------%
%------------------- BEGIN Function launchrv2oe.m ------------------------%
%-------------------------------------------------------------------------%
function oe = launchrv2oe(rv,vv,mu);

K  = [0;0;1];
hv = cross(rv,vv);
nv = cross(K,hv);
n  = sqrt(nv.'*nv);
h2 = (hv.'*hv);
v2 = (vv.'*vv);
r  = sqrt(rv.'*rv);
ev = 1/mu *( (v2-mu/r)*rv - (rv.'*vv)*vv );
p  = h2/mu;
e  = sqrt(ev.'*ev);
a  = p/(1-e*e);
i  = acos(hv(3)/sqrt(h2));
Om = acos(nv(1)/n);
if (nv(2)<0-eps),
  Om = 2*pi-Om;
end;
om = acos(nv.'*ev/n/e);
if (ev(3)<0),
  om = 2*pi-om;
end;
nu = acos(ev.'*rv/e/r);
if (rv.'*vv<0),
  nu = 2*pi-nu;
end;
oe = [a; e; i; Om; om; nu];

%-------------------------------------------------------------------------%
%---------------------- END Function launchrv2oe.m -----------------------%
%-------------------------------------------------------------------------%

