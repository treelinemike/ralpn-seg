% restart
close all; clear; clc;

% base directory for DICOM images
basePath = 'H:\CT\31584-008\CT 125797-125799 axial';

% text file with voxel coordinates of segmentation mask 
segFile = 'H:\CT\31584-008\CT 125797-125799 axial\RK_grayvalues.txt';

% extract list of DICOM files (all files w/o extensions)
% ref: https://www.mathworks.com/matlabcentral/answers/431023-list-all-and-only-files-with-no-extension
allFilesInDir = dir(basePath);
allFilenames = {allFilesInDir.name};
filesWithExtensionsMask = contains(allFilenames,'.');
allFilenames(filesWithExtensionsMask) = [];

% data storage
fileData = [];
allSegData = [];

% determine z location of each slice
for fileIdx = 1:length(allFilenames)
   
    thisFileFullPath = [basePath '\' allFilenames{fileIdx}];
    dinf = dicominfo(thisFileFullPath);
    fileData(fileIdx,:) = [dinf.ImagePositionPatient(3)];
end
[fileData,sortOrder] = sortrows(fileData,1);
fileData = [sortOrder fileData];   % [ fileIndex, DICOM slice #, Actual Z position in mm ]

% determine slice thickness
% this can come from any of the DICOM files, so we just use the last one
sliceThk = dinf.SliceThickness;  % in mm

% determine pixel spacing
pixSpace = dinf.PixelSpacing;
if(pixSpace(1) ~= pixSpace(2))
    error('Nonsquare pixels!');
end
pixSpace = pixSpace(1);

% load segmentation data
fid = fopen(segFile);
segData = textscan(fid,'%f%f%f%f','Delimiter',',','CollectOutput',1);
segData = segData{1};
fclose(fid);

% determine which slices we need
segData(:,3) = segData(:,3)-sliceThk/2;
sliceLocs = unique(segData(:,3),'stable');

% extract DICOM images and pair with the segmentations
for sliceIdx = 1:length(sliceLocs)
   
   % determine which DICOM image needs to be opened
   sliceToFind = sliceLocs(sliceIdx);
   thisFileIdxIdx = find(fileData(:,2) == sliceToFind);
   if( isempty(thisFileIdxIdx) )
       error('Slice not found!');
   end
   thisFileIdx = fileData(thisFileIdxIdx,1);
   thisFileFullPath = [basePath '\' allFilenames{thisFileIdx}];
   
   % load and store DICOM image
   img = dicomread(thisFileFullPath);
   minVal = -199;
   maxVal = 239;
%    minVal = double(min(img(:)));
%    maxVal = double(max(img(:)));
   img8 = uint8((img-minVal)*(255/(maxVal-minVal)));   
   allSegData(sliceIdx).img = img8;
   
   % store segmentation list
   seglist = segData(segData(:,3) == sliceToFind,1:2);
   allSegData(sliceIdx).seglist = seglist;
   seglistIdx = round(seglist/pixSpace)+repmat([257 256],size(seglist,1),1); % not sure why the x direction is off by a pixel, but 257 aligns MATLAB output to view in Mimics
    
   % compute segmentation mask labels
   imgLabels = uint8(zeros(size(img8)));
   for segPointIdx = 1:size(seglistIdx,1)
      imgLabels(seglistIdx(segPointIdx,2),seglistIdx(segPointIdx,1)) = 1;  % don't forget x is cols, y is rows!
   end
   allSegData(sliceIdx).labels = imgLabels;
   
   % generate masked image
   maskedImg = uint8(zeros(size(img8,1),size(img8,2),3));
   for layerIdx = 1:3
       maskedImg(:,:,layerIdx) = img8.*uint8(~imgLabels);
   end
   maskedImg(:,:,1) = maskedImg(:,:,1) + 255*imgLabels;
   allSegData(sliceIdx).maskedImg = maskedImg;
   
end

figure;
for sliceIdx = 1:length(allSegData)
    imshow( allSegData(sliceIdx).maskedImg );
    axis equal;
    title(sprintf('Slice Index %d @ %8.2f mm',sliceIdx,sliceLocs(sliceIdx)));
    drawnow;
    pause(0.1);
end
