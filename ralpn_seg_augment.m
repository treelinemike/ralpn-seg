function [im_out,seg_mask_out] = ralpn_seg_augment( im_in, seg_mask_in, ...
    contrastRange, ...
    angRange, ...
    transRange, ...
    warpCtrRange, ...
    warpSDMin, ...
    warpSDMax, ...
    maxWarpMag)

% set image and segmentation mask to input values
im = im_in;
seg_mask = seg_mask_in;

% contrast perturbation
% don't need to update segmentation mask for this
contrastFactor = contrastRange*rand(1)+1-0.5*contrastRange;
im = uint8(((double(im)/255).^contrastFactor)*255);

% warping perturbation
warpCtr = warpCtrRange*(rand(1,2)-0.5);
warpSD = (warpSDMax-warpSDMin)*rand(1)+warpSDMin;
x_range = (-256:255)/256;
y_range = (-256:255)/256;
[X,Y] = ndgrid(x_range,y_range);
Z = exp(-((X-warpCtr(1)).^2+(Y-warpCtr(2)).^2)/sqrt(warpSD));
[DX,DY] = gradient(Z,x_range,y_range);
mags = sqrt(DX.^2+DY.^2);
max_mag = max(mags(:));
DX = 1*round((maxWarpMag/max_mag)*DX);
DY = 1*round((maxWarpMag/max_mag)*DY); 
D = zeros(size(DX,1),size(DX,2),2);
D(:,:,1) = DX;
D(:,:,2) = DY;
im = imwarp(im,D);
seg_mask = imwarp(seg_mask,D);

% visualize warping
% figure;
% subplot(1,2,1);
% surf(X,Y,Z,'EdgeColor','none');
% view([0 90]);
% axis equal;
% xlim([-1 1]);
% ylim([-1 1]);
% subplot(1,2,2);
% quiver(X,Y,DX,DY);
% axis equal;
% xlim([-1 1]);
% ylim([-1 1]);

% rotational perturbation
thisAng = angRange*(rand(1)-0.5);   % [deg]
im = imrotate(im,thisAng,'nearest','crop');
seg_mask = imrotate(seg_mask,thisAng,'nearest','crop');

% translational perturbation
thisTrans = transRange*(rand(1,2)-0.5); % [px]
im = imtranslate(im,thisTrans,'nearest');
seg_mask = imtranslate(seg_mask,thisTrans,'nearest');

% collect results
im_out = im;
seg_mask_out = seg_mask;

end