% test of data augmentation function
% for RALPN U-Net based segmentation

% restart
close all; clear all; clc;

% define paths
image_path = 'C:\Users\f002r5k\GitHub\ralpn-seg\image\006_037_image.png';
seg_mask_path = 'C:\Users\f002r5k\GitHub\ralpn-seg\label\006_037_label.png';

% define segment colors
seg_colors = [ ...
    0.00 0.00 0.00; ... % background class
    1.00 0.75 0.00; ... % LK
    0.00 1.00 0.00; ... % RK
    1.00 0.33 0.33; ... % AA
    0.33 0.33 1.00; ... % IVC
    ];

% set input image and mask
im_in = imread(image_path);
seg_mask_in = imread(seg_mask_path);

% create random augmentation
tic;
[im_out, seg_mask_out] = ralpn_seg_augment( im_in, seg_mask_in, ...
    0.8, ...   % constrast range
    10, ...    % angle range [deg]
    30, ...    % translation range [pixels]
    4, ...     % warp center range in normalized image units (image dims -1 to 1)
    0.5, ...   % warp SD minimum 
    2.5, ...   % warp SD maximum
    20);       % max warp magnitude [pixels]
toc

% apply segmentation masks
im_in_masked = maskImage(rgb2gray(im_in),rgb2gray(seg_mask_in),seg_colors);
im_out_masked = maskImage(rgb2gray(im_out),rgb2gray(seg_mask_out),seg_colors);
im_masked_bad = maskImage(rgb2gray(im_out),rgb2gray(seg_mask_in),seg_colors);

% show results
figure;
set(gcf,'Position',[-4.600000e+00 2.074000e+02 1.407200e+03 0420]);
subplot(1,3,1);
imshow(im_in_masked);
title('\bfInput Image and Mask');

subplot(1,3,2);
imshow(im_masked_bad);
title('\bfOutput Image with Input Mask');

subplot(1,3,3);
imshow(im_out_masked);
title('\bfOutput Image and Mask');