function Y = Enhance_Polyfit31(Data,X)

n=3;
d=1;

[a, b]=size(X);
length=max([a b]);

x1=Data(1);
x2=Data(2);
x3=Data(3);

y(1)=Data(4);
y(2)=Data(5);
y(3)=Data(6);
y(4)=Data(7);
y(5)=Data(8);
y(6)=Data(9);

M=[[1,x1,x1^2,x1^3,x1^4,x1^5,];
[1,x2,x2^2,x2^3,x2^4,x2^5,];
[1,x3,x3^2,x3^3,x3^4,x3^5,];
[0,1,2*x1,3*x1^2,4*x1^3,5*x1^4,];
[0,1,2*x2,3*x2^2,4*x2^3,5*x2^4,];
[0,1,2*x3,3*x3^2,4*x3^3,5*x3^4,];
];

A=inv(M)*transpose(y);

Y=zeros(length,1);

for i=1:n*(d+1)
	Y(:)=Y(:)+A(i)*X(:).^(i-1);
end

return

