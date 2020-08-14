% Test conversion of Mimics CT segmentation to two 3D arrays to serve as input to a CNN (U-Net?):
%   * imageVolume: 512x512x128
%   * labelVolume: 512x512x128
%
% Once CT is segmented, in Mimics, right click on mask (or just in list of
% masks) and choose "Export Grayvalues...". Export the files in HU (doesn't
% really matter, but helpful for debugging), and save using conventions
% shown below in segFiles cell array.
%
% Also provide startSliceLoc and endSliceLoc which are the z coordinates of
% the center of the slice containing the inferior aspect of the L5 vertebrae
% (start), and the superior aspect of the T11 vertabrae (end)
%
% Script also writes a segmentation *.nrrd file which can be loaded (along
% with the corresponding DICOM stack) into 3D Slicer 4.11.0 or newer. Note
% that this requires NRRD/NHDR R/W from MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/66645-nrrd-nhdr-reader-and-writer
%
% Author:  M. Kokko
% Created: 28-Jun-2020

% restart
close all; clear; clc;

% options
doMakeVideo = 0;
doAnimate = 0;
doReslice = 0;
doSaveNRRD = 0;
doUNetExtract = 0;
doShowUNetImages = 0;
minHUScaleVal = -208;
maxHUScaleVal = 218;

% location of datasets along with start and end Z positions
% as defined by inferoior aspect of L5 (start), and superior aspect of T11
% (end); Cite Fananapazir2019 for justification of this range
dataSets = {
%     'H:\CT\31584-001',-2600,-1165; % 31584-001_full file!!
%     'H:\CT\31584-001',-1415,-1180; % 31584-001
    'H:\CT\31584-002\31584-003 6511 6514 CT',-388.69,-177.44;
    'H:\CT\31584-003\31584-003 6315 CT',-265.76,-54.35; % 31584-003
    'H:\CT\31584-004\31584-004-CT',513.2,753.2; % 31584-004
    'H:\CT\31584-005\CT 892882',1809.5,2028.5; % 31584-005
    'H:\CT\31584-006\CT 5761-5765',-269.4,-26.5; % 31584-006
    'H:\CT\31584-007\CT 9871-9873',-516.3,-267.3; % 31584-007
    'H:\CT\31584-008\CT 125797-125799 axial',1656.00,1896.00; % 31584-008
    'H:\CT\31584-009\CT 3061-3064',-1203.5,-973.5; % 31584-009
    'H:\CT\31584-010\CT 9001-9004',1435.00,1660.00; % 31584-010
    };

% text file with voxel coordinates of segmentation mask
segFiles = {...
    'LK_grayvalues.txt', ... % LK
    'RK_grayvalues.txt', ... % RK
    'AA_grayvalues.txt', ... % AA
    'IVC_grayvalues.txt', ... % IVC
    };

% define colors to use in masking...
segColors = [ ...
    0.00 0.00 0.00; ... % background class
    1.00 0.75 0.00; ... % LK
    0.00 1.00 0.00; ... % RK
    1.00 0.33 0.33; ... % AA
    0.33 0.33 1.00; ... % IVC
    ];

% storage for number of pixels in each class
classCounts = zeros(size(segColors,1),1);

% which file should we use right now?
% dataIdx = 2;
for dataIdx = 1:length(dataSets)
basePath = dataSets{dataIdx,1};
startSliceLoc = dataSets{dataIdx,2};
endSliceLoc = dataSets{dataIdx,3};

% extract list of DICOM files (all files w/o extensions)
% ref: https://www.mathworks.com/matlabcentral/answers/431023-list-all-and-only-files-with-no-extension
allFilesInDir = dir(basePath);
allFilenames = {allFilesInDir.name};
filesWithExtensionsMask = contains(allFilenames,'.');
allFilenames(filesWithExtensionsMask) = [];

%% determine z location of each slice
fileData = [];
for fileIdx = 1:length(allFilenames)
    thisFileFullPath = [basePath '\' allFilenames{fileIdx}];
    dinf = dicominfo(thisFileFullPath);
    fileData(fileIdx,:) = [dinf.ImagePositionPatient(3)];  % dinf.ImagePositionPatient(3) and dinf.SliceLocation should be identical!
end
[fileData,sortOrder] = sortrows(fileData,1);
fileData = [sortOrder fileData];   % [ file index in allFilenames, actual Z position in mm ]

% crop to specified ROI
sliceMask = (fileData(:,2) >= startSliceLoc) & (fileData(:,2) <= endSliceLoc);
fileData(~sliceMask,:) = [];

%% extract necessary image parameters
% these can come from any of the DICOM files, so we'll just use the last
% one that we opened
% first, determine slice SPACING (we don't want thickness here)
sliceSpacing = mean(diff(fileData(:,2)));
% sliceThk = dinf.SliceThickness;  % in mm

% now determine pixel spacing
pixSpace = dinf.PixelSpacing;
if(pixSpace(1) ~= pixSpace(2))
    error('Nonsquare pixels!');
end
pixSpace = pixSpace(1);

% finally, coordinates of the upper left pixel
% TODO: WE ASSUME THAT THIS IS CONSTANT FOR ENTIRE STACK!
% subtract off half the pixel width to get the coordinates of the UL
% corner
% need to add one pixel to x direction for consistency with MIMICS, don't
% entirely understand this yet?
ulPixCoords = dinf.ImagePositionPatient(1:2)'-pixSpace/2 + [-1 0]*pixSpace;

%% load and adjust all segmentation data
% segmentation z locations are taken on superior aspect of voxel
% whereas slice z locations are taken at center of voxel
% adjust all to be at center of voxel
maskData = [];
for segFileIdx = 1:length(segFiles)
    
    fid = fopen([basePath '\' segFiles{segFileIdx}]);
    thisSegData = textscan(fid,'%f%f%f%f','Delimiter',',','CollectOutput',1);
    thisSegData = thisSegData{1};
    fclose(fid);
    
    % determine which slices we need
    maskData(segFileIdx).slices = unique(thisSegData(:,3),'stable');
    
    % save all data in structure
    % [ voxelCtrZLocation, PixIdxX, PixIdxY ]
    pixelLocs =  round((thisSegData(:,1:2) - repmat(ulPixCoords,size(thisSegData,1),1))/pixSpace);
    maskData(segFileIdx).data = [thisSegData(:,3)-sliceSpacing/2 pixelLocs];
end

%% extract DICOM images and pair with the segmentations
allSegData = [];
for sliceIdx = 1:length(fileData)
    
    % get z location of this slice
    thisZLoc = fileData(sliceIdx,2);
    
    % determine which DICOM image needs to be opened
    thisFileIdx = fileData(sliceIdx,1);
    thisFileFullPath = [basePath '\' allFilenames{thisFileIdx}];
    
    % load and store DICOM image and first rescale into HU
    dinf = dicominfo(thisFileFullPath);
    img = double(dicomread(thisFileFullPath))*dinf.RescaleSlope + dinf.RescaleIntercept;  % will need coversion to uint8 to be standardized!
    
    % rescale HU to grayscale based on a "pretty good" mapping identified in
    % Mimics
    img8 = uint8((img-minHUScaleVal)*(255/(maxHUScaleVal-minHUScaleVal)));
    allSegData(sliceIdx).img = img8;
    
    % initialize segmentation mask
    seg_mask = zeros(size(img8));
    
    % prepare segmentation mask
    img8_masked = uint8(zeros(size(img8,1),size(img8,2),3));
    for layerIdx = 1:3
        img8_masked(:,:,layerIdx) = img8;
    end
    img8_masked_hsv = rgb2hsv(img8_masked);
    
    % generate overall segmentation mask for this slice
    for maskIdx = 1:length(maskData)
        
        % get coordinates of pixels to mask
        pixelLocs = maskData(maskIdx).data(  abs(maskData(maskIdx).data(:,1) - thisZLoc) < 0.1 ,2:3);
        
        % generate a mask for this image and this classification label
        thisMask = zeros(size(seg_mask));
        for pixelIdx = 1:size(pixelLocs,1)
            thisMask( pixelLocs(pixelIdx,2), pixelLocs(pixelIdx,1)) = 1;
        end
        
        % apply mask to overall mask
        seg_mask(thisMask ~= 0) = maskIdx;
        
    end
    
    % store segmentation mask and the masked image
    allSegData(sliceIdx).img8 = img8;
    allSegData(sliceIdx).seg_mask = seg_mask;
    
    % compute the masked image
    allSegData(sliceIdx).img8_masked = maskImage(img8,seg_mask,segColors);
    
    % store z location of this slice
    allSegData(sliceIdx).z_loc = dinf.ImagePositionPatient(3);
    
    % update counts
    for classIdx = 1:length(classCounts)
        classCounts(classIdx) = classCounts(classIdx) + nnz( seg_mask == (classIdx -1) );
    end
end

%% generate and export segmentation
segData = uint8(zeros(512,512, length(allSegData) ));
for i = 1:length(allSegData)
    segData(:,:,i) = allSegData(i).seg_mask;  % NOT TRANSPOSED YET, DO THIS LATER
end
delta_x = dinf.PixelSpacing(1);
delta_y = dinf.PixelSpacing(2);
delta_z = sliceSpacing;
spaceDir = [delta_x delta_y delta_z];
spaceDirMat = diag(spaceDir);


% write to file
if(doSaveNRRD)
% note: requires Slicer 4.11.0 or newer (loading nrrd crashes for Slicer
% 4.10.0, see:
% https://discourse.slicer.org/t/segment-editor-crashes-on-loaded-segments/9294/3)
headerInfo_new.content = 'matlab_export';
headerInfo_new.data = permute(segData,[2 1 3]);  % note permute does transpose! see: https://www.mathworks.com/matlabcentral/answers/162418-3-d-matrix-transpose
headerInfo_new.type = 'uint8';
headerInfo_new.dimension = 3;
headerInfo_new.space = 'left-posterior-superior';
headerInfo_new.sizes = size(segData);
for i = 1:3 
    headerInfo_new.spacedirections{i} = sprintf('(%19.17f,%19.17f,%19.17f)',spaceDirMat(:,i));
end
headerInfo_new.spacedirections_matrix = spaceDirMat;
headerInfo_new.kinds = {'domain'  'domain'  'domain'};
headerInfo_new.endian = 'little';
headerInfo_new.encoding = 'gzip';
headerInfo_new.spaceorigin = dinf.ImagePositionPatient;  % TODO UPDATE Z POSITION FOR STARTING SLICE!
headerInfo_new.spaceorigin(3) = fileData(1,2);
exportFilename = sprintf('seg_export_%03d.nrrd',dataIdx);
nhdr_nrrd_write(exportFilename, headerInfo_new, true);
disp(['Wrote ' sprintf('seg_export_%03d.nrrd',dataIdx)]);
end

%% now reslice to get a 512x512x128 voxel volume
% this might not be the best approach...
% ideally actual CT would be resliced and manually labeled
% but we'll try a simple reslicing in MATLAB first...
if(doReslice)
% make list of upper Z coordinate bounds for each slice
sliceUpperBounds = [allSegData.z_loc]'+pixSpace/2;

% generate new slice centers
newSliceCenters = linspace(allSegData(1).z_loc,allSegData(end).z_loc,128)';

% initialize data storage
imageVolume = zeros(512,512,128);
labelVolume = zeros(512,512,128);
allSegDataResampled = [];

% map image and label mask from original slicing
% note: this does NOT interpolate!
for newSliceIdx = 1:length(newSliceCenters)
    newSliceLoc = newSliceCenters(newSliceIdx);
    oldSliceToUse = find(sliceUpperBounds >= newSliceLoc,1,'first');
    thisImage = allSegData(oldSliceToUse).img8;
    thisMask = allSegData(oldSliceToUse).seg_mask;
    thisMaskedImg = allSegData(oldSliceToUse).img8_masked;
    
    % add to 3D data structures
    imageVolume(:,:,newSliceIdx) = thisImage;
    labelVolume(:,:,newSliceIdx) = thisMask;
    
    % add to a MATLAB struct for display
    allSegDataResampled(newSliceIdx).img8 = thisImage;
    allSegDataResampled(newSliceIdx).img8_masked = thisMaskedImg ;
    allSegDataResampled(newSliceIdx).seg_mask = thisMask;
    allSegDataResampled(newSliceIdx).z_loc = newSliceLoc;
end
end

%% produce animation and save video if desired
if(doAnimate)
% choose data to show
dispData = allSegData;
% dispData = allSegDataResampled;

% show each slice
figure;
for sliceIdx = 1:length(dispData)
    imshow( dispData(sliceIdx).img8_masked );
    axis equal;
    title(sprintf('Labeled Slice @ z = %8.2f mm',dispData(sliceIdx).z_loc));
    drawnow;
    
    if(doMakeVideo)
        thisImgFile = sprintf('frame%03d.png',sliceIdx);
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    else
        pause(0.1);
    end
end

% save animation as movie file
if(doMakeVideo)
    system(['ffmpeg -y -r 10 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 output.mp4']);
    system('del frame*.png');
end
end

%% extract some slices for testing U-Net
if(doUNetExtract)
% also show raw image, classification mask, and masked image
% 31584-001: 14-27
% 31584-002: 95-193
% 31584-003: 130-222
% 31584-004: 19-25
% 31584-005: 11-34
% 31584-006: 140-235
% 31584-007: 21-41
% 31584-008: 20-29
% 
startFrame = 1;
endFrame = length(allSegData);
% startFrame = 14;
% endFrame = 27;
numFrames = (endFrame-startFrame)+1;
ralpnData2D.image = uint8(zeros(512,512,numFrames));
ralpnData2D.label = uint8(zeros(512,512,numFrames));
if(doShowUNetImages)
    figure;
    set(gcf,'Position',[0169 0204 1375 0460]);
end
frameCount = 1;
for sliceIdx = startFrame:endFrame
    
    %     thisImage = imageVolume(:,:,sliceIdx);
    %     thisMask = labelVolume(:,:,sliceIdx);
    thisImage = allSegData(sliceIdx).img8;
    thisMask = allSegData(sliceIdx).seg_mask;
    
    ralpnData2D.image(:,:,frameCount) = thisImage;
    ralpnData2D.label(:,:,frameCount) = thisMask;
    frameCount = frameCount + 1;
    
    if(doShowUNetImages)
        img8_masked = maskImage(thisImage,thisMask,segColors);
        
        subplot(1,3,1);
        imshow(thisImage,[]);
        title('\bfRaw Image');
        
        subplot(1,3,2);
        imshow(thisMask+1,segColors);
        title('\bfClassification Mask');
        
        subplot(1,3,3);
        imshow(img8_masked);
        title('\bfMasked Image');
        drawnow;
    end

end
save(sprintf('ralpnData2D_%03d.mat',dataIdx),'ralpnData2D');
end

end

% compute class weights
% classCounts = [366699177; 1515535; 1623015; 613105; 745072]; % from cases 1-10
% from cases 2-10 (saving 1 for final check)
w = (1./classCounts)/sum(1./classCounts)
