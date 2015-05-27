function xyz_roipeak = roipeak(A,xyz,roi)
sz3d = size(A);
if(numel(size(roi)) == numel(sz3d) && all(size(roi)==sz3d))
    dmask = roi>0;
else
    dmask = roimask(sz3d,xyz,roi(1));
end
A(~dmask(:)) = nan;
[~,maxidx] = nanmax(A(:));
[xp, yp, zp] = ind2sub(sz3d,maxidx);
xyz_roipeak = [xp yp zp];