function L2 = rot_L2(ang)

L2 = [cos(ang) 0 -sin(ang); ...
             0 1         0; ...
      sin(ang) 0  cos(ang)];

return

