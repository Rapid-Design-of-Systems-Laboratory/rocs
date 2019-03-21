function out = xplt(NP,pop,vec,flag)
% xplt plots a coloured point with coordinates vec(1), vec(2)
% on a 3D-surface if flag is not equal 0. Otherwise NP colored
% points, stored in pop, will be plotted.
%
% Example: xplt(NP,pop)
%          where pop is a two-dimensional array of NP points
% 
%
% Used by: der.m

if (flag == 0)         %---draw entire population----------
  for i=1:NP
    x1=pop(i,1);
    x2=pop(i,2);
    z1=100*(x2-x1.^2).^2+(1-x1).^2;
    plot3(x1,x2,z1,'r.', ...
	'EraseMode','none', ...
	'MarkerSize',15);
    drawnow; %---Draws current graph now
    out = [];
  end
else
  x1 = vec(1);
  x2 = vec(2);
  z1=100*(x2-x1.^2).^2+(1-x1).^2;
  plot3(x1,x2,z1,'r.', ...
      'EraseMode','none', ...
      'MarkerSize',15);
  drawnow; %---Draws current graph now
  out = [];
end
