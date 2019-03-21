function Y = Enhance_Polyfit21(Data,X)

n=2;
d=1;

[a, b]=size(X);
length=max([a b]);

x1=Data(1);
x2=Data(2);

y(1)=Data(3);
y(2)=Data(4);
y(3)=Data(5);
y(4)=Data(6);

M=[[1,x1,x1^2,x1^3,];
[1,x2,x2^2,x2^3,];
[0,1,2*x1,3*x1^2,];
[0,1,2*x2,3*x2^2,];
];

A=inv(M)*transpose(y);

Y=zeros(length,1);

for i=1:n*(d+1)
	Y(:)=Y(:)+A(i)*X(:).^(i-1);
end

return

