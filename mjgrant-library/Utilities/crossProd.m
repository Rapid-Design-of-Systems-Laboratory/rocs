function c=crossProd(a,b)
%cross-product function (c=a x b)
%a,b,c are assumed 3-dimensional
c=[a(2)*b(3)-b(2)*a(3)
  -a(1)*b(3)+b(1)*a(3)
   a(1)*b(2)-b(1)*a(2)];
 
