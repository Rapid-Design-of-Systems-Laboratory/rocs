function createOptimizationPlots

close all; clc;

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

linewidth = 2;

x = linspace(-10,10,101);
y = linspace(-10,10,101);

[X,Y] = meshgrid(x,y);

J = X.^2 + Y.^2;

contourVals = linspace(min(min(J)),max(max(J)),11);

xConst = linspace(-10,0,101);
m = -1;
b = -9;
yConst = m*xConst+b;
JConst = xConst.^2 + yConst.^2;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Cost Function %%
%%%%%%%%%%%%%%%%%%%%%%%%

figure;
surf(X,Y,J,'EdgeColor','none','LineStyle','none');
hold on;
[c,h] = contour(X,Y,J,contourVals);
% clabel(c,h);
xlabel('X');
ylabel('Y');
zlabel('Z = X^2 + Y^2');
presentation_plot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Cost Function with Constraint %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
surf(X,Y,J,'EdgeColor','none','LineStyle','none');
hold on;
[c,h] = contour(X,Y,J,contourVals);
% clabel(c,h);
xlabel('X');
ylabel('Y');
zlabel('Z = X^2 + Y^2');
presentation_plot;
% alpha(0.4)

plot(xConst,yConst,'k','Linewidth',linewidth);
plot3(xConst,yConst,JConst,'k','Linewidth',linewidth);

return

