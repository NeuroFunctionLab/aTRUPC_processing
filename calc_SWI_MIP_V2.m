function [Im_MIP, Im_SWI4MIP_resliced, ind_SWI4MIP_surround_aTRUPC, Im_minIP] = calc_SWI_MIP(varargin)
% Input:
% dicom_path_aTRUPC: dicom folder of aTRUPC flow-comp images (the first series)
% dicom_path_SWI: dicom folder of SWI
% bCoreg: boolean flag, if bCoreg=1, do inplane coregistration between SWI and aTRUPC
% Im_aTRUPC: flow-comp aTRUPC image reconstructed from rawdata (i.e., mag1)
% Output:
% Im_MIP: maximum intensity projection of SWI at the same spatial location as aTRUPC
% Im_SWI4MIP_resliced: resliced SWI images used to calculate Im_MIP
% ind_SWI4MIP_surround_aTRUPC: slice index of the two Im_SWI4MIP_resliced images that "sandwich" the location of aTRUPC 
% Im_minIP: minimal intensity projection showing the veins
% Dengrong Jiang, 9/26/2023
dicom_path_aTRUPC = varargin{1};
dicom_path_SWI = varargin{2};
if nargin > 2
    bCoreg = varargin{3};
    Im_aTRUPC = varargin{4};
else
    bCoreg = 0;
end
FileList_aTRUPC = dir([dicom_path_aTRUPC, filesep, '*.IMA']);
Filename_list_aTRUPC = extractfield(FileList_aTRUPC,'name');
aTRUPC_Info = dicominfo(strcat(dicom_path_aTRUPC, filesep, Filename_list_aTRUPC{1}));
aTRUPC_ImagePosition = aTRUPC_Info.ImagePositionPatient;
aTRUPC_ImageOrientation = aTRUPC_Info.ImageOrientationPatient;
slcthick_aTRUPC = aTRUPC_Info.SliceThickness;
pixel_spacing_aTRUPC = aTRUPC_Info.PixelSpacing; % [row_spacing, column_spacing]
delta_j_aTRUPC = pixel_spacing_aTRUPC(1);
delta_i_aTRUPC = pixel_spacing_aTRUPC(2);

[Im_SWI, SWI_ImageOrientation, SWI_ImagePosition_array, delta_i_SWI, delta_j_SWI, slcthick_SWI] = read_3D_dicom(dicom_path_SWI);
[nRow_SWI, nCol_SWI, nSlc_SWI] = size(Im_SWI);
dist_SWI_aTRUPC = sqrt(sum((SWI_ImagePosition_array - repmat(aTRUPC_ImagePosition, [1, nSlc_SWI])).^2, 1));
[~, sorted_ind_dist_SWI] = sort(dist_SWI_aTRUPC, 'ascend');
% aTRUPC actually has the same position as one of the SWI mIP images, each
% mIP is generated from 8 SWI images. The position of the mIP image is the
% average position of the 4th and 5th SWI images used to calculate the mIP
% Therefore, the position of aTRUPC is actually in the middle of two
% SWI images
ind_SWI_surround_aTRUPC = sorted_ind_dist_SWI(1:2);
ind_SWI_surround_aTRUPC = sort(ind_SWI_surround_aTRUPC, 'ascend');
SWI_mIP_ImagePosition = mean(SWI_ImagePosition_array(:, ind_SWI_surround_aTRUPC), 2);

num_SWI4MIP = ceil(slcthick_aTRUPC/slcthick_SWI);
ind_SWI4MIP_surround_aTRUPC = [round(num_SWI4MIP/2), round(num_SWI4MIP/2)+1];
Im_SWI4MIP = Im_SWI(:, :, ind_SWI_surround_aTRUPC(1)-round(num_SWI4MIP/2)+1 : ind_SWI_surround_aTRUPC(1)-round(num_SWI4MIP/2)+num_SWI4MIP);

%% reslice based on DICOM tags
cntIm = dicomread(strcat(dicom_path_aTRUPC, filesep, Filename_list_aTRUPC{1}));
[nRow_aTRUPC, nCol_aTRUPC] = size(cntIm);
% Compute the matrices to convert [i;j;1] to [x;y;1] for both source and target
% Note that for SWI, the average position to the two SWI images that
% "sandwich" the aTRUPC image is used, in other words, use the SWI mIP
% image that has the same location as aTRUPC
Matrix_source = [SWI_ImageOrientation(1)*delta_i_SWI, SWI_ImageOrientation(4)*delta_j_SWI, SWI_mIP_ImagePosition(1);
                 SWI_ImageOrientation(2)*delta_i_SWI, SWI_ImageOrientation(5)*delta_j_SWI, SWI_mIP_ImagePosition(2);
                 0, 0, 1;
                 ];

Matrix_target = [aTRUPC_ImageOrientation(1)*delta_i_aTRUPC, aTRUPC_ImageOrientation(4)*delta_j_aTRUPC, aTRUPC_ImagePosition(1);
                 aTRUPC_ImageOrientation(2)*delta_i_aTRUPC, aTRUPC_ImageOrientation(5)*delta_j_aTRUPC, aTRUPC_ImagePosition(2);
                 0, 0, 1;
                 ];
% Compute the matrix to convert target [i_t; j_t; 1] to source [i_s; j_s; 1]
Matrix_target_to_source = Matrix_source\Matrix_target;

% reslice
[x_source, y_source] = meshgrid([0:1:nCol_SWI-1], [0:1:nRow_SWI-1]);

[x_target, y_target] = meshgrid([0:1:nCol_aTRUPC-1], [0:1:nRow_aTRUPC-1]);

coordinate_in_source = Matrix_target_to_source*[x_target(:).';...
                                           y_target(:).';...
                                           ones(1, numel(x_target))];
x_target_in_source = reshape(coordinate_in_source(1, :), [nRow_aTRUPC, nCol_aTRUPC]);
y_target_in_source = reshape(coordinate_in_source(2, :), [nRow_aTRUPC, nCol_aTRUPC]);

Im_SWI4MIP_resliced = zeros(nRow_aTRUPC, nCol_aTRUPC, num_SWI4MIP);
for iSlc = 1:num_SWI4MIP
    Image_source = Im_SWI4MIP(:,:,iSlc);
    Im_SWI4MIP_resliced(:,:,iSlc) = interp2(x_source, y_source, Image_source,...
                     x_target_in_source, y_target_in_source, 'linear', 0);
end

%% inplane coregistration
if bCoreg
    Ref2d = imref2d([nRow_aTRUPC, nCol_aTRUPC]);
    Ref3d = imref3d([nRow_aTRUPC, nCol_aTRUPC, num_SWI4MIP]);
    
    [optimizer_MI,metric_MI] = imregconfig("multimodal");
    optimizer_MI.InitialRadius = 0.001;
    optimizer_MI.Epsilon = 1.5e-4;
    optimizer_MI.GrowthFactor = 1.01;
    optimizer_MI.MaximumIterations = 300;
    metric_MI.NumberOfHistogramBins = 25;
    
%     the position of aTRUPC is actually in the middle of two SWI images
    Im_moving = mean(Im_SWI4MIP_resliced(:,:,ind_SWI4MIP_surround_aTRUPC), 3);
    Im_SWI4MIP_resliced_ori = Im_SWI4MIP_resliced;
    
    tform2D = imregtform(Im_moving,Ref2d,Im_aTRUPC,Ref2d,"rigid",optimizer_MI,metric_MI);
    
    T_2D = tform2D.T;
    tMat = [T_2D(1:2, 1:2), zeros(2);
            0, 0, 1, 0;
            T_2D(3, 1:2), 0, 1];
    tform = affine3d(tMat);
    Im_SWI4MIP_resliced = imwarp(Im_SWI4MIP_resliced_ori, tform, 'OutputView', Ref3d);
end

Im_MIP = max(Im_SWI4MIP_resliced, [], 3);
Im_minIP = min(Im_SWI4MIP_resliced, [], 3);
