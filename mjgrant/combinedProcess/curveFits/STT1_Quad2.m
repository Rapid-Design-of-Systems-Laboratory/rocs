function Y = STT1_Quad2(Data,X)

n=2;

[a, b]=size(X);
length=max([a b]);

x1=Data(1);
x2=Data(2);

y1=Data(3);
y2=Data(4);
y3=Data(5);

for i=1:length
X1=X(i);
Y(i)=y1 + y2*(x2/2 - x1/2 + (X1 - x2)^2/(2*(x1 - x2))) + (y3*(x1 - x2))/2 + y2*(x1/2 - x2/2) - (y3*(X1 - x1)^2)/(2*(x1 - x2));
end

return

