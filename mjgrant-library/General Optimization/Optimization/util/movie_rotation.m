function movie_rotation

axis vis3d;
aviobj = avifile('pareto_rotation_movie.avi','Compression','Indeo5','FPS',8);

for counter = 1 : 1 : 90
  camorbit(1,0);
  frame = getframe(gcf);
  aviobj = addframe(aviobj,frame); % Add figure to movie
end

aviobj = close(aviobj);

return

