function mask = roimask(sz3d,xyz,roisize)
[my, mx, mz] = meshgrid(1:sz3d(2),1:sz3d(1),1:sz3d(3));
d = sqrt((mx-xyz(1)).^2+(my-xyz(2)).^2+(mz-xyz(3)).^2);
mask = d<=roisize;
