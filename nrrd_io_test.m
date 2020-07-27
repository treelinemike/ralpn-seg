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
    subplot(1,2,1);
    cla
%     imshow(thisFrame',[]);  % NOTE THE TRANSPOSE!!!
    imshow(maskImage(ralpnData2D.image(:,:,frameIdx),thisFrame',segColorsSlicer));  % NOTE THE TRANSPOSE!!!

    title(['\bfFrame ' num2str(frameIdx) ' - Slicer']);
    
    subplot(1,2,2);
    cla
    imshow(maskImage(ralpnData2D.image(:,:,frameIdx),ralpnData2D.label(:,:,frameIdx),segColorsMimics));  % NOTE THE TRANSPOSE!!!
    title(['\bfFrame ' num2str(frameIdx) ' - Mimics']);
    
    drawnow;
end