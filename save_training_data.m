% restart
close all; clear; clc;

% get files to load
basePath = '.';
fileList = 	cellstr(ls([basePath '/' 'ralpnData2D_*.mat']));

mkdir('./image/');
mkdir('./label/');

for fileIdx = 1:length(fileList)
    thisFileName = fileList{fileIdx};
    tok = regexp(thisFileName,'ralpnData2D_(\d+).mat','tokens');
    ctIdx = str2double(tok{1}{1});
    load([basePath '/' thisFileName]);
    numExamples = size(ralpnData2D.image,3);
    for exampIdx = 1:numExamples
        image = repmat(ralpnData2D.image(:,:,exampIdx),1,1,3);
        label = repmat(ralpnData2D.label(:,:,exampIdx),1,1,3);
        imageFile = sprintf('%03d_%03d_image.png',ctIdx,exampIdx);
        labelFile = sprintf('%03d_%03d_label.png',ctIdx,exampIdx);
        imwrite(uint8(image),[basePath '/image/' imageFile]);
        imwrite(uint8(label),[basePath '/label/' labelFile]);
    end
end