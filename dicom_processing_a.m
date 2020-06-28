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
% Author:  M. Kokko
% Created: 28-Jun-2020

% restart
close all; clear; clc;

% options
doMakeVideo = 0;

% base directory for DICOM images
% basePath = 'H:\CT\31584-008\CT 125797-125799 axial';
basePath = 'G:\CT\31584-007\CT 9871-9873';

% slice locations to start and end
% startSliceLoc = 1656.00; % inferior aspect of on L5
% endSliceLoc = 1896.00; % superior aspect of T11
startSliceLoc = -516.3; % inferior aspect of on L5 (cite Fananapazir2019 for justification of L5 to T11 range)
endSliceLoc = -267.3; % superior aspect of T11

% text file with voxel coordinates of segmentation mask
segFiles = {...
    'LK_grayvalues.txt', ... % LK
    'RK_grayvalues.txt', ... % RK
    'AA_grayvalues.txt', ... % AA
    'IVC_grayvalues.txt', ... % IVC
    };

% define colors to use in masking...
segColors = [ ...
    0.00 1.00 0.00; ... % LK
    0.00 1.00 0.00; ... % RK
    1.00 0.33 0.33; ... % AA
    0.33 0.33 1.00; ... % IVC
    ];

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
    fileData(fileIdx,:) = [dinf.ImagePositionPatient(3)];
end
[fileData,sortOrder] = sortrows(fileData,1);
fileData = [sortOrder fileData];   % [ file index in allFilenames, actual Z position in mm ]

% crop to specified ROI
sliceMask = (fileData(:,2) >= startSliceLoc) & (fileData(:,2) <= endSliceLoc);
fileData(~sliceMask,:) = [];

%% extract necessary image parameters
% these can come from any of the DICOM files, so we'll just use the last
% one that we opened
% first, determine slice thickness
sliceThk = dinf.SliceThickness;  % in mm

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
ulPixCoords = dinf.ImagePositionPatient(1:2)'-pixSpace/2 + [-1 0];

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
    maskData(segFileIdx).data = [thisSegData(:,3)-sliceThk/2 pixelLocs];
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
    minVal = -208;
    maxVal = 218;
    %    minVal = double(min(img(:)));
    %    maxVal = double(max(img(:)));
    img8 = uint8((img-minVal)*(255/(maxVal-minVal)));
    allSegData(sliceIdx).img = img8;
    
    % initialize segmentation mask
    seg_mask = zeros(size(img8));
    
    % prepare masked image
    img8_masked = uint8(zeros(size(img8,1),size(img8,2),3));
    for layerIdx = 1:3
        img8_masked(:,:,layerIdx) = img8;
    end
    img8_masked_hsv = rgb2hsv(img8_masked);
    
    % generate overall segmentation mask for this slice
    for maskIdx = 1:length(maskData)
        
        % get coordinates of pixels to mask
        pixelLocs = maskData(maskIdx).data(maskData(maskIdx).data(:,1) == thisZLoc,2:3);
        
        % generate a mask for this image and this classification label
        thisMask = zeros(size(seg_mask));
        for pixelIdx = 1:size(pixelLocs,1)
            thisMask( pixelLocs(pixelIdx,2), pixelLocs(pixelIdx,1)) = 1;
        end
        
        % apply mask to overall mask
        seg_mask(thisMask ~= 0) = maskIdx;
        
        % apply mask to masked image
        thisColorHSV = rgb2hsv(segColors(maskIdx,:));
    
        % update hue and saturation
        img8_masked_hsv(:,:,1) = img8_masked_hsv(:,:,1) + thisColorHSV(1)*thisMask;
        img8_masked_hsv(:,:,2) = img8_masked_hsv(:,:,2) + thisColorHSV(2)*thisMask;
        
    end
    
    % store segmentation mask and the masked image
    allSegData(sliceIdx).img8 = img8;
    allSegData(sliceIdx).seg_mask = seg_mask;
    allSegData(sliceIdx).img8_masked = hsv2rgb(img8_masked_hsv);
    
    % store z location of this slice
    allSegData(sliceIdx).z_loc = dinf.ImagePositionPatient(3);
    
end

%% now reslice to get a 512x512x128 voxel volume
% this might not be the best approach...
% ideally actual CT would be resliced and manually labeled
% but we'll try a simple reslicing in MATLAB first...

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

%% produce animation and save video if desired

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