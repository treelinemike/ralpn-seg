% restart
close all; clear all; clc;

image_path = 'C:\Users\f002r5k\GitHub\ralpn-seg\image\001_007_image.png';
seg_mask_path = 'C:\Users\f002r5k\GitHub\ralpn-seg\label\001_007_label.png';

seg_colors = [ ...
    0.00 0.00 0.00; ... % background class
    1.00 0.75 0.00; ... % LK
    0.00 1.00 0.00; ... % RK
    1.00 0.33 0.33; ... % AA
    0.33 0.33 1.00; ... % IVC
    ];

im = imread(image_path);
seg_mask = imread(seg_mask_path);
im_masked = maskImage(rgb2gray(im),rgb2gray(seg_mask),seg_colors);

% adjust contrast
contrastRange = 0.8;
contrastFactor = contrastRange*rand(1)+1-0.5*contrastRange
im = uint8(((double(im)/255).^contrastFactor)*255);

thisAng = -5+10*rand(1)   % [deg]
im = imrotate(im,thisAng,'nearest','crop');
seg_mask = imrotate(seg_mask,thisAng,'nearest','crop');

thisTrans = -20+40*rand(1,2) % [px]
im = imtranslate(im,thisTrans,'nearest');
seg_mask = imtranslate(seg_mask,thisTrans,'nearest');


SD = 2*rand(1)+0.5
ctr = -2+4*rand(1,2)

x_range = (-256:255)/256;
y_range = (-256:255)/256;

[X,Y] = ndgrid(x_range,y_range);
Z = exp(-((X-ctr(1)).^2+(Y-ctr(2)).^2)/sqrt(SD));
[DX,DY] = gradient(Z,x_range,y_range);
mags = sqrt(DX.^2+DY.^2);
max_mag = max(mags(:));
set_max_mag = 20;
DX = 1*round((set_max_mag/max_mag)*DX);
DY = 1*round((set_max_mag/max_mag)*DY); 
D = zeros(size(DX,1),size(DX,2),2);
D(:,:,1) = DX;
D(:,:,2) = DY;
im_warped = imwarp(im,D);
seg_mask_warped = imwarp(seg_mask,D);
im_masked_warped = maskImage(rgb2gray(im_warped),rgb2gray(seg_mask_warped),seg_colors);
im_masked_warped_bad = maskImage(rgb2gray(im_warped),rgb2gray(seg_mask),seg_colors);

figure;
subplot(1,2,1);
surf(X,Y,Z,'EdgeColor','none');
subplot(1,2,2);
quiver(X,Y,DX,DY);
axis equal;
figure;
set(gcf,'Position',[-4.600000e+00 2.074000e+02 1.407200e+03 0420]);
subplot(1,3,1);
imshow(im_masked);
subplot(1,3,2);
imshow(im_masked_warped_bad);
subplot(1,3,3);
imshow(im_masked_warped);



