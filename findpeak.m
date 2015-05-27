function xyz_peak = findpeak(A, xyz0)


try
    
    %error('die');
%smooth image for better gradient descent
Asmooth = smooth3(A,'gaussian', [5 5 5]);
Asmooth(isnan(Asmooth)) = 0;
valfunc = @(p)(-interpn(Asmooth,p(1),p(2),p(3),'linear'));

options = optimset('Algorithm','interior-point','display','off');

x_range = [10 10 10];
xyz_peak = [];

%if starting voxel is zero-value, find nearest non-zero to start
if(valfunc(xyz0) == 0)
    sz3d = size(A);
    [my, mx, mz] = meshgrid(1:sz3d(2),1:sz3d(1),1:sz3d(3));
    d2 = (mx-xyz0(1)).^2+(my-xyz0(2)).^2+(mz-xyz0(3)).^2;
    nzidx = find(Asmooth(:) > 0);
    [~,minidx] = min(d2(nzidx));
    xyz0 = [mx(nzidx(minidx)) my(nzidx(minidx)) mz(nzidx(minidx))];
end


%iterate gradient descent several times to make sure we've found a local peak
% (since we only let it search a limited range around starting point)
fprintf('Finding flip angle peak');
for iter = 1:5
    fprintf('...');    

    xmax = max(min([xyz0+x_range],size(A)),[1 1 1]);
    xmin = max(min([xyz0-x_range],size(A)),[1 1 1]);
    xyz_peak = fmincon(valfunc,xyz0,[],[],[],[],xmin,xmax,[],options);

    %stop once our result is the same as our starting point
    if(all(xyz0 == xyz_peak))
        break;
    else
        xyz0 = xyz_peak;
    end

end

fprintf(' Done\n');

catch err
    if(strcmpi(err.identifier,'MATLAB:license:checkouterror'))
        fprintf('\nOptimization Toolbox license unavailable. Unable to run findpeak.\n');
    end
	xyz_peak = xyz0;
end
