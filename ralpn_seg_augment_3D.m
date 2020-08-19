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

% convert one hot to dense
[~,seg_mask_dense]=max(seg_mask,[],4);
seg_mask_dense = seg_mask_dense-1;

% rotational perturbation
thisAng = angRange*(rand(1)-0.5);   % [deg]
im = imrotate3(im,thisAng,[0 0 1],'nearest','crop');
seg_mask_dense = imrotate3(seg_mask_dense,thisAng,[0 0 1],'nearest','crop');

% warping perturbation
warpCtr = warpCtrRange*(rand(1,3)-0.5);
warpSD = (warpSDMax-warpSDMin)*rand(1)+warpSDMin;
x_range = (-1*ceil(size(im,1)/2):(size(im,1)-ceil(size(im,1)/2)-1))/ceil(size(im,1));
y_range = (-1*ceil(size(im,2)/2):(size(im,2)-ceil(size(im,2)/2)-1))/ceil(size(im,2));
z_range = (-1*ceil(size(im,3)/2):(size(im,3)-ceil(size(im,3)/2)-1))/ceil(size(im,3));
[X,Y,Z] = ndgrid(x_range,y_range,z_range);
F = exp(-0.5*((X-warpCtr(1)).^2+(Y-warpCtr(2)).^2+(Z-warpCtr(3)).^2)/sqrt(warpSD));
[DX,DY,DZ] = gradient(F,x_range,y_range,z_range);
mags = sqrt(DX.^2+DY.^2+DZ.^2);
max_mag = max(mags(:));
DX = round((maxWarpMag/max_mag)*DX);
DY = round((maxWarpMag/max_mag)*DY); 
DZ = round((maxWarpMag/max_mag)*DZ); 
D = zeros(size(DX,1),size(DX,2),size(DX,3),3);
D(:,:,:,1) = DX;
D(:,:,:,2) = DY;
D(:,:,:,3) = DZ;
im = imwarp(im,D,'FillValues',0);
seg_mask_dense = imwarp(seg_mask_dense,D,'nearest');

% translational perturbation
thisTrans = transRange*(rand(1,3)-0.5); % [px]
thisTrans(3) = 0; % don't shift in z!
im = imtranslate(im,thisTrans,'nearest');
seg_mask_dense = imtranslate(seg_mask_dense,thisTrans,'nearest');

% covert dense to one hot
seg_mask = zeros(size(seg_mask));
for labelIdx = 1:size(seg_mask,4)
    seg_mask(:,:,:,labelIdx) = (seg_mask_dense == labelIdx-1);
end

% collect results
im_out = im;
seg_mask_out = seg_mask;

end