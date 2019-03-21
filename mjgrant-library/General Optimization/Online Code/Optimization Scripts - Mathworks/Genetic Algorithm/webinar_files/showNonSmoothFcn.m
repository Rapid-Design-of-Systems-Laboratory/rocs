function showNonSmoothFcn(fcn,range)
if(nargin == 0)
    fcn = @rastriginsfcn;
    range = [-5,5;-5,5];
end

pts = 25;
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

surf(x,y,values)
shading interp
light
lighting phong
hold on
contour(x,y,values)
rotate3d
view(37,60)

%Annotations
figure1 = gcf;
% Create arrow
annotation1 = annotation(figure1,'arrow',[0.5946 0.4196],[0.9024 0.6738]);
% Create textbox
annotation2 = annotation(...
  figure1,'textbox',...
  'Position',[0.575 0.9071 0.1571 0.07402],...
  'FitHeightToText','off',...
  'FontWeight','bold',...
  'String',{'Start point'});
 
% Create textarrow
annotation3 = annotation(...
  figure1,'textarrow',...
  [0.3679 0.4661],[0.1476 0.3214],...
  'String',{'Non-differentiable regions'},...
  'FontWeight','bold');
 
% Create arrow
annotation4 = annotation(figure1,'arrow',[0.1196 0.04107],[0.1381 0.5429]);
 
% Create textarrow
annotation5 = annotation(...
  figure1,'textarrow',...
  [0.7411 0.5321],[0.05476 0.1381],...
  'LineWidth',2,...
  'Color',[1 0 0],...
  'String',{'Smooth region'},...
  'FontWeight','bold',...
  'TextLineWidth',2,...
  'TextEdgeColor',[1 0 0]);
 
% Create arrow
annotation6 = annotation(...
  figure1,'arrow',...
  [0.8946 0.9179],[0.05714 0.531],...
  'Color',[1 0 0]);
 
