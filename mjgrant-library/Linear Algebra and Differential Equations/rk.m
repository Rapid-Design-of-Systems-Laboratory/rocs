function [X, Y] = rk(x, y, xf, n)

h = (xf - x) / n;
X = x;
Y = y;

for i = 1:n
    k1 = f(x,y);
    k2 = f(x+h/2,y+h*k1/2);
    k3 = f(x+h/2, y+h*k2/2);
    k4 = f(x+h, y+h*k3);
    y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    x = x + h;
    X = [X;x];
    Y = [Y;y];
end