% restart
close all; clear; clc;

% get files to load
basePath = '.';
fileList = 	cellstr(ls([basePath '/' 'ralpnData2D_*.mat']));

mkdir('./image/');
mkdir('./seg_mask/');

for fileIdx = 1:length(fileList)
    thisFileName = fileList{fileIdx};
    tok = regexp(thisFileName,'ralpnData2D_(\d+).mat','tokens');
    ctIdx = str2double(tok{1}{1});
    load([basePath '/' thisFileName]);
    numExamples = size(ralpnData2D.image,3);
    for exampIdx = 1:numExamples
        
        % extract image from .mat file
        image = repmat(ralpnData2D.image(:,:,exampIdx),1,1,3);
        
%         % look at interior of image
%         image2 = uint8(zeros(size(image)));
%         for layerIdx = 1:3
%            image2(97:416,97:416,layerIdx) =  image(97:416,97:416,layerIdx);
%         end
%         imshow(image2);
%         drawnow;

        % extract FULL segmentation mask from .mat file
        seg_mask_full = repmat(ralpnData2D.label(:,:,exampIdx),1,1,3);

        % produce cropped segmentation mask
        seg_mask_crop = zeros(320,320,3);
        for layerIdx = 1:3
           seg_mask_crop(:,:,layerIdx) =  seg_mask_full(97:416,97:416,layerIdx);
        end
        
        % decide which segmentation mask to use
        seg_mask = seg_mask_crop;

        % write images to file
        imageFile = sprintf('%03d_%03d_image.png',ctIdx,exampIdx);
        segMaskFile = sprintf('%03d_%03d_seg_mask.png',ctIdx,exampIdx);
        imwrite(uint8(image),[basePath '/image/' imageFile]);
        imwrite(uint8(seg_mask),[basePath '/seg_mask/' segMaskFile]);
        
        
    end
end