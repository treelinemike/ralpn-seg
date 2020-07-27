% restart
close all; clear; clc

% colormaps
segColorsMimics = [ ...
    0.00 0.00 0.00; ... % background class
    1.00 0.75 0.00; ... % LK
    0.00 1.00 0.00; ... % RK
    1.00 0.33 0.33; ... % AA
    0.33 0.33 1.00; ... % IVC
    ];
segColorsSlicer = [ ...
    0.00 0.00 0.00; ... % background class
    1.00 0.75 0.00; ... % LK
    0.00 1.00 0.00; ... % RK
    ];
segColorsDiff = [ ...
    0.00 0.00 0.00; ... % background class
    0.33 0.33 1.00; ... % Mimics Only
    1.00 0.33 0.33; ... % Slicer Only
    ];

% read SLICER data
headerInfo = nhdr_nrrd_read('C:\Users\f002r5k\Desktop\ralpn_local\slicer-test\kidney_segs_a.nrrd', true);

% read MIMICS data
load('ralpnData2D_001_all.mat');

% check size of data
if( size(headerInfo.data,3) ~= size(ralpnData2D.image,3))
    error('Size mismatch!');
end

figure;

for frameIdx = 1:size(headerInfo.data,3)
    thisFrame = headerInfo.data(:,:,frameIdx);
    subplot(1,3,1);
    cla
    %     imshow(thisFrame',[]);  % NOTE THE TRANSPOSE!!!
    slicerMask = thisFrame';
    imshow(maskImage(ralpnData2D.image(:,:,frameIdx),slicerMask,segColorsSlicer));  % NOTE THE TRANSPOSE!!!

    title(['\bfFrame ' num2str(frameIdx) ' - Slicer']);
    
    subplot(1,3,2);
    cla
    mimicsMask = ralpnData2D.label(:,:,frameIdx);
    imshow(maskImage(ralpnData2D.image(:,:,frameIdx),mimicsMask,segColorsMimics));  % NOTE THE TRANSPOSE!!!
    title(['\bfFrame ' num2str(frameIdx) ' - Mimics']);

    subplot(1,3,3);
    cla
    diffMask = double(mimicsMask > 0 & mimicsMask < 3) - double(slicerMask > 0);
    diffMask( diffMask == -1 ) = 2;
    diffMask = uint8(diffMask);
    imshow(maskImage(ralpnData2D.image(:,:,frameIdx),diffMask,segColorsDiff));  % NOTE THE TRANSPOSE!!!
    title(['\bfFrame ' num2str(frameIdx) ' - Mask Diff']);

    
    drawnow;
end