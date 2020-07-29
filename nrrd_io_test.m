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
set(gcf,'Position',[0089 2.194000e+02 1.404800e+03 4.064000e+02]);
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

dicomFilename = 'C:\Users\f002r5k\Desktop\ralpn_local\slicer-test\31584-001-ct\9043';
dicomHeader = dicominfo(dicomFilename);
delta_x = dicomHeader.PixelSpacing(1);
delta_y = dicomHeader.PixelSpacing(2);
delta_z = 5; %TODO: fix this! use mean(diff(dicomHeader.SliceLocation))
spaceDir = [delta_x delta_y delta_z];
spaceDirMat = diag(spaceDir);


% write to file
% note: requires Slicer 4.11.0 or newer (loading nrrd crashes for Slicer
% 4.10.0, see:
% https://discourse.slicer.org/t/segment-editor-crashes-on-loaded-segments/9294/3)
headerInfo_new.content = 'matlab_export';
headerInfo_new.data = permute(ralpnData2D.label,[2 1 3]);  % note permute does transpose! see: https://www.mathworks.com/matlabcentral/answers/162418-3-d-matrix-transpose
headerInfo_new.type = 'uint8';
headerInfo_new.dimension = 3;
headerInfo_new.space = 'left-posterior-superior';
headerInfo_new.sizes = size(ralpnData2D.label);
for i = 1:3 
    headerInfo_new.spacedirections{i} = sprintf('(%19.17f,%19.17f,%19.17f)',spaceDirMat(:,i));
end
headerInfo_new.spacedirections_matrix = spaceDirMat;
headerInfo_new.kinds = {'domain'  'domain'  'domain'};
headerInfo_new.endian = 'little';
headerInfo_new.encoding = 'gzip';
headerInfo_new.spaceorigin = dicomHeader.ImagePositionPatient;
nhdr_nrrd_write('testMatlabA.nrrd', headerInfo_new, true);