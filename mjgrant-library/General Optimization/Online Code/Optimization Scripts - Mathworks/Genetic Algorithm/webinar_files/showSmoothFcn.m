function showSmoothFcn(fcn,range)


pts = 100;
span = diff(range')/(pts - 1);
x = range(1,1): span(1) : range(1,2);
y = range(2,1): span(2) : range(2,2);

pop = zeros(pts * pts,2);
k = 1;
for i = 1:pts
    for j = 1:pts
        pop(k,:) = [x(i),y(j)];
        k = k + 1;
    end
end

values = feval(fcn,pop);
values = reshape(values,pts,pts);
clf;
surf(x,y,values)
shading interp

light
lighting phong
hold on
rotate3d
view(37,60)
set(gcf,'Renderer','opengl');
set(gca,'ZLimMode','manual');