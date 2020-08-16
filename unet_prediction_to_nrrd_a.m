% CONVERT TENSORFLOW OUTPUT TO NRRD
% FOR IMPORT INTO 3D SLICER

% restart
close all; clear; clc;

% load already-processed masks, etc. for this example
load('ralpnData2D_001.mat');

% load predicted masks - this is output from Tensorflow / U-Net
load('final_masks.mat');
segData = uint8(final_masks);  % need to conver to uint8

% EXPORT NRRD FILE!
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
    headerInfo_new.spacedirections{i} = sprintf('(%19.17f,%19.17f,%19.17f)',ralpnData2D.spaceDirMat(:,i));
end
headerInfo_new.spacedirections_matrix = ralpnData2D.spaceDirMat;
headerInfo_new.kinds = {'domain'  'domain'  'domain'};
headerInfo_new.endian = 'little';
headerInfo_new.encoding = 'gzip';
headerInfo_new.spaceorigin = ralpnData2D.spaceorigin;
exportFilename = 'auto_seg.nrrd';
nhdr_nrrd_write(exportFilename, headerInfo_new, true);
disp(['Wrote ' exportFilename]);