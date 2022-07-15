% Convert MIMICS segmentations into a format importable into 3D Slicer
% Author: Mike Kokko
% Updated: 15-Jul-2022

% NOTES:
% * DICOM folder must contain **ONLY** DICOM images
% * Once CT is segmented, in Mimics, right click on mask (or just in list of
%   masks) and choose "Export Grayvalues...". Export the files in HU (doesn't
%   really matter, but helpful for debugging), and save using conventions
%   shown below in segFiles cell array.
% * Script writes a segmentation *.nrrd file which can be loaded (along
%   with the corresponding DICOM stack) into 3D Slicer 4.11.0 or newer. Note
%   that this requires NRRD/NHDR R/W from MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/66645-nrrd-nhdr-reader-and-writer

% restart
close all; clear all; clc;

% options
exportFilename = 'manual_seg.nrrd';

dcm_dir = 'G:\.shortcut-targets-by-id\1DXab6fRFPLxd5fj49tY3mc-NZnjbf0IM\STUDY Lymph Node\CTs\A1\slicer_conversion_test\dcm';
gv_filenames = {'R Final_grayvalues_mak.txt'};
gv_colors = [   1 0 0  ];

% dcm_dir = 'C:\Users\f002r5k\Desktop\ct\dcm_dir';
% gv_filenames = {'C:\Users\f002r5k\Desktop\ct\other\AA_grayvalues.txt','D:\CT\31584-001\IVC_grayvalues.txt'};
% gv_colors = [   1 0 0; ...
%                 0 0 1;  ];

% load grayvalues files into tables
gv = [];
disp('Loading grayvalues tables...');
for gv_idx = 1:length(gv_filenames)
    gv(gv_idx).tab = readtable(gv_filenames{gv_idx});
end

% load all DICOM data (ref: https://www.mathworks.com/matlabcentral/answers/431023-list-all-and-only-files-with-no-extension?s_tid=mwa_osa_a)
dcm_data = [];
slice_locs = [];
disp('Loading DICOM data...');
files = dir(fullfile(dcm_dir, '*'));
filenames = {files.name};
filenames(ismember(filenames, {'.', '..'})) = [];
for file_idx = 1:length(filenames)
    disp(num2str(file_idx));
    filename_full = fullfile(dcm_dir,filenames{file_idx});
    dcm_data(file_idx).info = dicominfo(filename_full);
    dcm_data(file_idx).im = dicomread(filename_full);
    slice_locs(file_idx) =  dcm_data(file_idx).info.ImagePositionPatient(3);
end

% sort slice locations (b/c filenames may not be ordered correctly)
[slice_locs,slice_order] = sort(slice_locs);
slice_inc = mean(diff(slice_locs));  % why isn't this a DICOM field???

% now create masks
segdata = uint8(zeros(size(dcm_data(1).im,1),size(dcm_data(1).im,2),length(slice_order)));
for order_idx = 1:length(slice_order)

    % get transformation
    file_idx = slice_order(order_idx);
    z_ctr = dcm_data(file_idx).info.ImagePositionPatient(3);
    z_gv = z_ctr + 0.5*slice_inc;
    X = dcm_data(file_idx).info.ImageOrientationPatient(1:3);
    Y = dcm_data(file_idx).info.ImageOrientationPatient(4:6);
    delta_i = dcm_data(file_idx).info.PixelSpacing(2);  % COLUMN spacing
    delta_j = dcm_data(file_idx).info.PixelSpacing(1);  % ROW spacing
    S = dcm_data(file_idx).info.ImagePositionPatient;
    M = [X*delta_i, Y*delta_j, zeros(3,1), S; 0 0 0 1];
    Mxy = [M(1:2,1:2) M(1:2,4); 0 0 1]; % maps [row col 1]' to [x y 1]'

    % mask all segments
    slice_mask = zeros(size(dcm_data(file_idx).im));
    dcm_img8 = uint8(255*double(dcm_data(file_idx).im)/double(max(dcm_data(file_idx).im(:))));
    im_overlay = uint8(zeros(size(dcm_data(file_idx).im,1),size(dcm_data(file_idx).im,2),3));
    for i = 1:3
        im_overlay(:,:,i) = dcm_img8;
    end

    for seg_idx = 1:length(gv)
        gv_tabmask = (abs(gv(seg_idx).tab.Var3 - z_gv) < 0.1);
        gv_subtab = gv(seg_idx).tab(gv_tabmask,:);

        gv_subtab.Var1 = gv_subtab.Var1 - 0.5*delta_j;
        gv_subtab.Var2 = gv_subtab.Var2 + 0.5*delta_i;
        gv_subtab.Var3 = gv_subtab.Var3 - 0.5*slice_inc;

        gv_xy = [gv_subtab.Var1 gv_subtab.Var2 ones(size(gv_subtab,1),1)]';
        gv_ij = round(inv(Mxy)*gv_xy);
        gv_ij = gv_ij(1:2,:) + [1;1];   % [1;1] b/c MATLAB is one indexed!

        gv_immask = zeros(size(dcm_data(file_idx).im));
        mask_idx = sub2ind(size(dcm_data(file_idx).im),gv_ij(2,:),gv_ij(1,:));
        gv_immask(mask_idx) = 1;

        % add to segdata
        % note: does not allow overlapping segments!
        % last segments have highest priority!
        segdata(:,:,order_idx) = uint8(segdata(:,:,order_idx).*uint8(~gv_immask)) + uint8(gv_immask.*seg_idx);

        % update overlay for display only
        im_overlay = imoverlay(im_overlay,gv_immask,gv_colors(seg_idx,:));
    end
    
    imshow(im_overlay);
    drawnow;
end

%%
% note: requires Slicer 4.11.0 or newer (loading nrrd crashes for Slicer 4.10.0, see:
% https://discourse.slicer.org/t/segment-editor-crashes-on-loaded-segments/9294/3)
delta_x = dcm_data(1).info.PixelSpacing(1);
delta_y = dcm_data(1).info.PixelSpacing(2);
delta_z = slice_inc;

% spaceDir = [delta_x delta_y delta_z];
spaceDir = [M(1,1) M(2,2) slice_inc];
spaceDirMat = diag(spaceDir);
spaceorigin = dcm_data(1).info.ImagePositionPatient;
spaceorigin(3) = slice_locs(1);              % UPDATE Z POSITION FOR STARTING SLICE!

headerInfo_new.content = 'matlab_export';
headerInfo_new.data = permute(segdata,[2 1 3]);  % note permute does transpose! see: https://www.mathworks.com/matlabcentral/answers/162418-3-d-matrix-transpose
headerInfo_new.type = 'uint8';
headerInfo_new.dimension = 3;
headerInfo_new.space = 'left-posterior-superior';  % TODO: SHOULD WE BE UPDATING THIS?!?
headerInfo_new.sizes = size(segdata);
for i = 1:3
    headerInfo_new.spacedirections{i} = sprintf('(%19.17f,%19.17f,%19.17f)',spaceDirMat(:,i));
end
headerInfo_new.spacedirections_matrix = spaceDirMat;
headerInfo_new.kinds = {'domain'  'domain'  'domain'};
headerInfo_new.endian = 'little';
headerInfo_new.encoding = 'gzip';
headerInfo_new.spaceorigin = spaceorigin;
nhdr_nrrd_write(exportFilename, headerInfo_new, true);
disp(['Wrote ' exportFilename]);
