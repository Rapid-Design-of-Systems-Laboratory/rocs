function L1 = rot_L1(ang)

L1 = [1         0        0; ...
      0  cos(ang) sin(ang); ...
      0 -sin(ang) cos(ang)];

return

