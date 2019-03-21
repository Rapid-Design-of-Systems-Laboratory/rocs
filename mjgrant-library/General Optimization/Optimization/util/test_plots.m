function test_plots

clear all; pack; close all; clc;

figure;
colordef(gcf,'black')
cla
load wind
[cx cy cz] = meshgrid(linspace(71,134,10),linspace(18,59,10),3:4:15);
daspect([1 1 1])
h=coneplot(x,y,z,u,v,w,cx,cy,cz,y,3);
set(h,'EdgeColor', 'none');
colormap(hsv);
box on;
axis tight
camproj perspective;
camva(35);
campos([175 10 85]);
camtarget([105 40 0])
camlight left;
lighting gouraud

return


data = load('results.mat');
in = data.in;
od = data.od;
clear data;

x_vals = 10;
y_vals = 10;
z_vals = 10;

x = od.fit_arch(1,:,od.iter);
y = od.fit_arch(2,:,od.iter);
z = od.fit_arch(3,:,od.iter);
v = od.pos_arch(1,:,od.iter);

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);
z_min = min(z);
z_max = max(z);

x_set = [x_min : (x_max-x_min)/x_vals : x_max];
y_set = [y_min : (y_max-y_min)/y_vals : y_max];
z_set = [z_min : (z_max-z_min)/z_vals : z_max];

% [xi,yi,zi] = meshgrid(x_set,y_set,z_set);
% w = griddata3(x,y,z,v,xi,yi,zi);

[xi,yi] = meshgrid(x_set,y_set);
w = griddata(x,y,z,xi,yi);

surf(xi,yi,w);

% p = patch(isosurface(xi,yi,zi,w,0.8));
% isonormals(xi,yi,zi,w,p);
% set(p,'FaceColor','blue','EdgeColor','none');
% view(3), axis equal, axis off, camlight, lighting phong

return

