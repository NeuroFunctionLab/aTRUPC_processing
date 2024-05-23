%% Note:
% 202306: change the ROI drawing of the static tissue, 
% 20230710: updated version
% 20230720: V9 : updated version, change the pixel indexing to physical
% distance for ROI drawing
% 20230801: V10: updates epsilon calculation for motion-contaminated repetitions
% 20240522: update Physical length option in adjustable parameters section
%%
close all;
clear;
tic

%% Project-specific Parameters
addpath('C:\Program Files\MATLAB\R2023b\toolbox\spm12'); % path of SPM12
eTE = [0.44, 40, 80]; % effective TEs, in ms
im_orientation_PE = 'sag'; % Fixed. 'sag' (sagittal) or 'ax_ap' (axial, PE direction=A>>P)

% ROI in sagittal view
NumROI = 6; % Fixed. Number of ROIs for quantifying T2 and Yv
NumPixInROI = [100, 100, 100, 40, 40, 40]; % number of pixels selected for T2 fitting in each ROI
Physcdist_1 = 30; % Unit: mm. Physical length of each ROI on SSS and SS. eg: 30mm for adults.
Physcdist_2 = 10; % Unit: mm. Physical length of the interval between the great vein and internal cerebral vein (starting from the ending point)

%% --------------- No change after this line ---------------

%% Load data
[filename, pathname] = uigetfile('*.dat', 'Select the raw data .dat file');
[~,NAME_file,fileEXT] = fileparts(filename);
hct = str2double(inputdlg('Hct'));

savefilename = 'results_real_imag.mat';
results_folder = [pathname, NAME_file, '_real_imag_semiAuto_V12'];
assert(~exist(results_folder, 'dir'), 'Folder exists!');
    
mkdir(results_folder);

spm_get_defaults;
global defaults;

%% Advanced Parameters
ref_img = []; % target image for realignment. When it's empty, use the 1st dynamic as the target
% use stored masks
bUseStoredROI = 0; % use stored ROI for assessing T2 and Yv
bUseStoredTissueROI = 0;
bLoadPreVesselMask = 0; % preliminary mask for vessel masking
% Yv fitting
bT2Yvmap = 1; % fit T2 and Yv maps
NumPixInROI_max = round(NumPixInROI*10); % max number of pixels used in each ROI
bT2Correction = 1;
T2fit_method = 'multiS0'; % avg, multiS0
vessel_thresh_ratio = 0.8; % thresholding to refine vessel mask
bUseMaxPixsInROI = 1; % if true, only use the max few pixels in the ROI
T2map_window_size = [9, 9];
ind2fit = [1:length(eTE)];
% ind2fit = [1:2];
% image reconstruction
addpath('GRAPPA_JDR');
bAdaptiveComb = 1; % if true, use adaptive coil combination when computing CD
bSquareImage = 1; % make a square image by zero-padding
bROfilter = 0; % apply hanning filter along kx
GRAPPA_kersize = [5,5];
GRAPPA_lambda = 0.01;
kspace_scale = 1e6;
deg_poly_eddycurrent = 2; % degree of polynomial for hyperplane fitting to correct for eddy current
% k-space Hanning filter
NslopePE = 25;
NslopeRO = 224;
% for display
MagAxisSame = [0,50];

%% T2 correction formula
Fun_T2Cor = @(x) 1.2002*x - 10.6276; % fitting by minimizing modified chi-square function, on 8 adults with 6 ROIs

%% Read raw data
mapVBVD_results = mapVBVD_DJ([pathname, filename], 'imaref', 1);
kspace_mapVBVD = mapVBVD_results{2}.image();
slice_quaternion = mapVBVD_results{2}.image.slicePos(:,1);
kspace_mapVBVD = single(kspace_mapVBVD); % to save memory, also because BART uses single (32bit) precision
kspace_mapVBVD = permute(kspace_mapVBVD, [1,3,4,2,10,6,9,5,7,8]);
kspace_mapVBVD = squeeze(kspace_mapVBVD);
[nRO, nLin, nCh, nSet, neTE, nRep] = size(kspace_mapVBVD);
kspace_header = mapVBVD_results{2}.hdr();
% inplane_rot_rad = kspace_header.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot;
ksp_param = kspace_header.MeasYaps.sKSpace;
MaxwellCoef = kspace_header.Meas.aflMaxwellCoefficients;
nPEy = ksp_param.lPhaseEncodingLines;
kspace_multi = single(zeros(nRO, nPEy, nCh, nSet, neTE, nRep));
kspace_multi(:, 1:nLin, :, :, :, :) = kspace_mapVBVD;
clear kspace_mapVBVD; % to save memory
voxelsize = [kspace_header.MeasYaps.sSliceArray.asSlice{1, 1}.dReadoutFOV/ksp_param.lBaseResolution,...
    kspace_header.MeasYaps.sSliceArray.asSlice{1, 1}.dPhaseFOV/(ksp_param.lBaseResolution*ksp_param.dPhaseResolution),...
    kspace_header.MeasYaps.sSliceArray.asSlice{1, 1}.dThickness];
VENC = kspace_header.MeasYaps.sAngio.sFlowArray.asElm{1, 1}.nVelocity;
Rfactor = kspace_header.MeasYaps.sPat.lAccelFactPE;
slice_center_pos = kspace_header.MeasYaps.sSliceArray.asSlice{1, 1}.sPosition;
slice_normal_vector = kspace_header.MeasYaps.sSliceArray.asSlice{1, 1}.sNormal;
% clear kspace_header; % to save memory
% clear mapVBVD_results; % to save memory
[nRO, nPE, nCh, nSet, neTE, nRep] = size(kspace_multi);
kspace_multi = kspace_multi.*kspace_scale;
if Rfactor > 1
    bGRAPPA = 1; % 1D GRAPPA in kx-ky plane
    RefKy = kspace_header.MeasYaps.sPat.lRefLinesPE;
else
    bGRAPPA = 0;
end

%% Remove incomplete repetition
cntK = kspace_multi(:,:,:,:,:,1);
numksample = sum(abs(cntK(:)) > 0);
for iRep = 2:nRep
    cntK = kspace_multi(:,:,:,:,:,iRep);
    if numksample*0.9 > sum(abs(cntK(:)) > 0) % detect incomplete repetition
        fprintf('Rep #%d is incomplete!\n', iRep);
        % remove incomplete repetition
        kspace_multi(:,:,:,:,:,iRep:end) = [];
        nRep = size(kspace_multi, 6);
        break;
    end
end

%% GRAPPA reconstruction and filtering
% k-space filter
if bROfilter
    NflatRO = nRO - 2*NslopeRO;
    [fltRO, fwhmRO] = HanningFilter1D(nRO, NflatRO, NslopeRO);
end

if bSquareImage % zero-padding to make an square image
    % when the actual number of phase encoding lines exceeds
    % base resolution * phase resolution, the Siemens sequence seems to
    % reduce delta-ky to keep the ky coverage as desired, and let the
    % FOVy to be increased, which can be cut-off after image is reconstructed
    % Thus, ky should be zero-padded to be larger than nRO/2 by the
    % factor of nPE/(BaseResolution*PhaseResolution)
    voxelsize = [voxelsize(1), voxelsize(1), voxelsize(3)];
    fprintf('Recon voxelsize = %.2fx%.2fx%.2fmm\n', voxelsize(1), voxelsize(2), voxelsize(3));
    nPE_square = floor( nRO/2 * nPE / (ksp_param.dPhaseResolution*ksp_param.lBaseResolution) );
    NflatPE = nPE-2*NslopePE;
    [fltPE, fwhmPE] = HanningFilter1D(nPE_square, NflatPE, NslopePE);
    fltPE = repmat(fltPE, [nRO, 1]);
    if bROfilter
        fltRO = repmat(fltRO.', [1, nPE_square]);
    end
end

Image_coil = cell(nRep, neTE, nSet);

for iRep = 1:nRep
    for i_eTE = 1:neTE
        for iSet = 1:nSet
            if bGRAPPA % Do GRAPPA recon one by one
                RefKx = RefKy;
                cntKspace = squeeze(kspace_multi(:,:,:,iSet,i_eTE,iRep));
                kCalib = crop(cntKspace, [RefKx, RefKy, nCh]);
                kspace_raw = GRAPPA_JDR_V2(cntKspace,kCalib,GRAPPA_kersize,GRAPPA_lambda, Rfactor);
            else
                kspace_raw = squeeze(kspace_multi(:,:,:,iSet,i_eTE,iRep));
            end
            
            kspace_square = zeros(nRO, nPE_square, nCh);
            kspace_square(:, floor((nPE_square+1)/2)-floor(nPE/2)+1:floor((nPE_square+1)/2)-floor(nPE/2)+nPE, :) = kspace_raw;
            
            for iCh = 1:nCh
                kspace_square(:,:,iCh) = kspace_square(:,:,iCh).*fltPE;
                if bROfilter
                    kspace_square(:,:,iCh) = kspace_square(:,:,iCh).*fltRO;
                end
            end
            
            im_coil = fftshift(fftshift(fft(ifft(ifftshift(ifftshift(kspace_square, 1),2),[],1),[],2),1),2);
            if bSquareImage % remove oversampling in RO and also the excessive FOV in PE
                im_square = im_coil(floor(nRO/4)+1:floor(nRO/4)+floor(nRO/2), floor((nPE_square+1)/2)-floor(nRO/4)+1:floor((nPE_square+1)/2)-floor(nRO/4)+floor(nRO/2), :, :);
            else
                % remove oversampling in RO
                im_square = im_coil(floor(nRO/4)+1:floor(nRO/4)+floor(nRO/2), :, :, :);
            end
            if strcmpi(im_orientation_PE, 'ax_ap')
                Image_coil{iRep, i_eTE, iSet} = permute(im_square, [2,1,3]);
%                 Image_coil{iRep, i_eTE, iSet} = flip(im_square, 1);
            elseif strcmpi(im_orientation_PE, 'sag')
                Image_coil{iRep, i_eTE, iSet} = im_square;
            end
        end
    end
end
[Nrow, Ncol, ~] = size(Image_coil{1});

%% Exclude motion-contaminated repetitions | 20230808
% save([results_folder, filesep, 'varforrep.mat']);
% epsilon = zeros(size(Image_coil,1),1);
% epsilon = exclude_rep_by_epsilon(results_folder);

% Exclude by fixed threshold for epsilon
% epi_thresh = 0.08;
% rep2exclude = find(epsilon>epi_thresh);
% if length(rep2exclude)==4
%     rep2exclude = input('Exclude repetitions?:');
%     Image_coil(rep2exclude, :, :) = []; % mark 
%     nRep = nRep - numel(rep2exclude);
% else
%     Image_coil(rep2exclude, :, :) = []; % mark 
%     nRep = nRep - numel(rep2exclude);
% end
Im2disp = cell(nRep, 1);
figure;
for iRep = 1:nRep
    subplot(2, nRep, iRep);
    Im2disp{iRep} = zeros(Nrow, Ncol);
    for iSet = 1:nSet - 1
        Im2disp{iRep} = Im2disp{iRep} + sum(abs(Image_coil{iRep,1,1} - Image_coil{iRep,1,iSet+1}).^2, 3);
    end
    Im2disp{iRep} = sqrt(Im2disp{iRep});
    imshow(Im2disp{iRep}, [0, 150]);
    title(sprintf('Rep #%d', iRep));
    subplot(2, nRep, iRep+nRep);
    imshow(sqrt(sum(abs(Image_coil{iRep,1,1}).^2, 3)), []);
    set(gca, 'fontsize', 15);
end

rep2exclude = input('Exclude repetitions?:');
Image_coil(rep2exclude, :, :) = []; % mark 
nRep = nRep - numel(rep2exclude);

handle_suplot=get(gcf,'children');
% The get(gcf, ‘children’) function returns the handles for the sub plots 
% in the reverse order (i.e., the last subplot as the first handle)
handle_suplot=handle_suplot(end:-1:1);
for iExclude = 1:length(rep2exclude)
    iRep = rep2exclude(iExclude);
    title(handle_suplot(iRep), sprintf('Rep #%d (Excluded)', iRep));
end
saveas(gcf, [results_folder, filesep, 'inspect_images.fig']);

%% Realign across repetitions
answer1 = questdlg(sprintf('Realign?'), ...
    'Realign?', ...
    'Yes','No','Yes');
if strcmp(answer1, 'Yes')
    Image_coil_ori = Image_coil; % save a copy of original images before registration
    % use eTE0 mag1 image of each repetition for co-registration
    im2reg = [];
    mag1_reg = cell(nRep, 1);
    for iRep = 1:nRep
        mag1_reg{iRep} = abs(CoilAdaptiveCombine_V2(Image_coil{iRep,1,1}, 1, 9, 0, eye(nCh)));
        if iRep < 10
            mag1_reg_file = [results_folder, filesep, sprintf('mag1_Rep0%d.img', iRep)];
        else
            mag1_reg_file = [results_folder, filesep, sprintf('mag1_Rep%d.img', iRep)];
        end
        write_hdrimg(mag1_reg{iRep}, mag1_reg_file, voxelsize, 16, 1);
        im2reg = [im2reg; mag1_reg_file];
    end
    
    if ~isempty(ref_img) % realign to a reference image
        target_img = [results_folder, filesep, 'mag1_Rep00'];
        [ref_path, ref_name, ref_ext] = fileparts(ref_img);
        copyfile([ref_path, filesep, ref_name, '.img'], [target_img, '.img']);
        copyfile([ref_path, filesep, ref_name, '.hdr'], [target_img, '.hdr']);
        im2reg = [[target_img, '.img']; im2reg];
    end
    
    % SPM realign: a least squares approach and a 6 parameter (rigid body) spatial transformation
    defs = defaults.realign;
    FlagsC = struct('quality',defs.estimate.quality,...
        'fwhm',5,'rtm',0);
    spm_realign(im2reg, FlagsC);
    which_writerealign = 2;
    mean_writerealign = 0;
    FlagsR = struct('interp',defs.write.interp,...
        'wrap',defs.write.wrap,...
        'mask',defs.write.mask,...
        'which',which_writerealign,'mean',mean_writerealign);    
    spm_reslice(im2reg, FlagsR);
    
    % read transformation parameters
    [target_path, target_name, ~] = fileparts(im2reg(1,:));
    motion_vec = importdata([target_path, filesep, 'rp_', target_name, '.txt']);
    
    % Note: Matlab rotation is around the origin (0,0), so to perform
    % rotation around the center of image, we need to first translate the
    % image center to the origin, do the rotation, then translate the image
    % back to its original center

    % to make the codes faster, here we treat each multi-channel
    % complex image as a 3D image, and generate a 3D Matlab affine
    % transform that only has x- and y-translation and z-rotation
    Ref3d = imref3d([Nrow, Ncol, nCh]);
    tX = mean(Ref3d.XWorldLimits);
    tY = mean(Ref3d.YWorldLimits);
    tZ = mean(Ref3d.ZWorldLimits);
    tTranslate2Origin = [1 0 0 0; 0 1 0 0; 0 0 1 0; -tX -tY -tZ 1]; 
    tTranslateBack  = [1 0 0 0; 0 1 0 0; 0 0 1 0; tX tY tZ 1];
    
    if ~isempty(ref_img) % realign to a reference image
        for iRep = 1:nRep
            % generate Matlab transform based on SPM transformation parameters
            % Note: the first column of the SPM motion vector is y-translation
            % and the sign is opposite to Matlab convention!
            % the second is x-translation

            % to make the codes faster, here we treat each multi-channel
            % complex image as a 3D image, and generate a 3D Matlab affine
            % transform that only has x- and y-translation and z-rotation
            theta = motion_vec(iRep+1, 6);
            tRot = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
             % Note: the first column of the SPM motion vector is y-translation
             % and the sign is opposite to Matlab convention!
             % the second is x-translation
            tYtrans = -motion_vec(iRep+1, 1)/voxelsize(1);
            tXtrans = motion_vec(iRep+1, 2)/voxelsize(2);
            tTrans = [1 0 0 0; 0 1 0 0; 0 0 1 0; tXtrans tYtrans 0 1];
            tMat = tTranslate2Origin*tRot*tTranslateBack*tTrans;
            tform = affine3d(tMat);

            % apply the transform to all images
            for i_eTE = 1:neTE
                for iSet = 1:nSet
                    im_real = imwarp(real(Image_coil_ori{iRep, i_eTE, iSet}), tform, 'OutputView', Ref3d);
                    im_imag = imwarp(imag(Image_coil_ori{iRep, i_eTE, iSet}), tform, 'OutputView', Ref3d);
                    Image_coil{iRep, i_eTE, iSet} = im_real + sqrt(-1)*im_imag;
                end
            end
        end
    else % realign to the first repetiion
        mag1_trans{1} = mag1_reg{1};

        for iRep = 2:nRep
            % generate Matlab transform based on SPM transformation parameters
            % Note: the first column of the SPM motion vector is y-translation
            % and the sign is opposite to Matlab convention!
            % the second is x-translation

            % to make the codes faster, here we treat each multi-channel
            % complex image as a 3D image, and generate a 3D Matlab affine
            % transform that only has x- and y-translation and z-rotation
            theta = motion_vec(iRep, 6);
            tRot = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
             % Note: the first column of the SPM motion vector is y-translation
             % and the sign is opposite to Matlab convention!
             % the second is x-translation
            tYtrans = -motion_vec(iRep, 1)/voxelsize(1);
            tXtrans = motion_vec(iRep, 2)/voxelsize(2);
            tTrans = [1 0 0 0; 0 1 0 0; 0 0 1 0; tXtrans tYtrans 0 1];
            tMat = tTranslate2Origin*tRot*tTranslateBack*tTrans;
            tform = affine3d(tMat);

            % apply the transform to all images
            for i_eTE = 1:neTE
                for iSet = 1:nSet
                    im_real = imwarp(real(Image_coil_ori{iRep, i_eTE, iSet}), tform, 'OutputView', Ref3d);
                    im_imag = imwarp(imag(Image_coil_ori{iRep, i_eTE, iSet}), tform, 'OutputView', Ref3d);
                    Image_coil{iRep, i_eTE, iSet} = im_real + sqrt(-1)*im_imag;
                end
            end
        end
    end
end

%% Complex average across repetitions
Image_coil_avg = cell(neTE, nSet);
for i_eTE = 1:neTE
    for iSet = 1:nSet
        Image_coil_avg{i_eTE, iSet} = zeros(size(Image_coil{1}));
    end
end

for i_eTE = 1:neTE
    for iSet = 1:nSet
        for iRep = 1:nRep
            Image_coil_avg{i_eTE, iSet} = Image_coil_avg{i_eTE, iSet} + Image_coil{iRep, i_eTE, iSet};
        end
        Image_coil_avg{i_eTE, iSet} = Image_coil_avg{i_eTE, iSet}/nRep;
    end
end

%% Adaptive combine multi-coil complex images
Image_combined = cell(neTE, nSet);
% use the eTE0 flow-compensated images to estimate the weights
[Image_combined{1,1}, weights, wtFull] = CoilAdaptiveCombine_V2(Image_coil_avg{1,1}, 1, 9, 1, eye(nCh));
% use the same weights for other images
for i_eTE = 1:neTE
    for iSet = 1:nSet
        if ~(i_eTE == 1 && iSet == 1)
            Image_combined{i_eTE, iSet} = sum(Image_coil_avg{i_eTE, iSet} .* conj(weights), 3);
        end
    end
end

%% Calculate complex difference for subsequent eddy current correction
CD = cell(neTE, nSet-1);
mag1 = cell(neTE, 1);
mag2 = cell(neTE, nSet-1);
dphi = cell(neTE, nSet-1);
Z1Z2product = cell(neTE, nSet-1);
for i_eTE = 1:neTE
    mag1{i_eTE} = abs(Image_combined{i_eTE, 1});
    for iSet = 1:nSet-1
        mag2{i_eTE, iSet} = abs(Image_combined{i_eTE, iSet+1});
        Z1Z2product{i_eTE, iSet} = Image_combined{i_eTE, iSet+1}.*conj(Image_combined{i_eTE, 1});
        dphi{i_eTE, iSet} = argument(Z1Z2product{i_eTE, iSet});
        CD{i_eTE, iSet} = abs(Image_combined{i_eTE, iSet+1} - Image_combined{i_eTE, 1});
    end
end

CD2mask = zeros(size(CD{1}));
for iSet = 1:nSet-1
    CD2mask = CD2mask + CD{1, iSet}.^2;
end
CD2mask = sqrt(CD2mask);

%% Segmentation & eddy current correction of static tissue mask | Auto-selection
% Use anatomical image:
AI = abs(Image_combined{1, 1});
BW_brain = extract_brain_mask(AI,0);
AI_2 = AI./max(AI(:)); 
C2 = labeloverlay(AI_2,BW_brain);
figure,imshow(C2);title('Brain ROI-Automatic');
saveas(gcf, [results_folder, filesep, 'pre_Auto_Brain_ROI.fig']);

% threshold to segment out the vessel
CD2mask_brain = CD2mask.*BW_brain; 
CD2mask_brain_2 = CD2mask_brain./max(CD2mask_brain(:)); 
BW_vessel = CD2mask_brain_2 > graythresh(CD2mask_brain_2(logical(BW_brain)));
CD_image = CD2mask./max(CD2mask(:));
C1 = labeloverlay(CD_image,BW_vessel);
figure,imshow(C1);title('pre-Static vessel ROI-Automatic');
saveas(gcf, [results_folder, filesep, 'pre_Auto_Vessel_ROI.fig']);

% substract to get static tissue image.
BW_static = BW_brain - BW_vessel;
% save the ROI
C = labeloverlay(CD_image,BW_static);
figure,imshow(C);title('pre-Static tissue ROI-Automatic');
saveas(gcf, [results_folder, filesep, 'pre_Auto_static_ROI.fig']);

%% Automatic eddy current correction
CD_correct = cell(neTE, nSet-1);
dphi_correct = cell(neTE, nSet-1);
phi_fit = cell(nSet-1, 1);
cc_fit = cell(nSet-1, 1);
Niter = 5;
bBatch = 1;
for iSet = 1:nSet-1
    [CD_correct(:, iSet), dphi_correct(:, iSet), phi_fit{iSet}, cc_fit{iSet}, mask_static_refined] = CDPhaCor_Iter_2D_V3...
        (mag1, mag2(:, iSet), dphi(:, iSet), Z1Z2product(:, iSet), BW_static, BW_brain, deg_poly_eddycurrent, Niter, bBatch, 1);
end

C3 = labeloverlay(CD_image,mask_static_refined);
figure,imshow(C3);title(['Static tissue ROI-Automatic, iter = ' num2str(Niter)]);
saveas(gcf, [results_folder, filesep, 'Refined_st_ROI_Auto.fig']);

CD_complex = cell(neTE, nSet-1);
for i_eTE = 1:neTE
    for iSet = 1:nSet-1
        CD_complex{i_eTE, iSet} = Image_combined{i_eTE, iSet+1}.*exp(-sqrt(-1).*phi_fit{iSet}) - Image_combined{i_eTE, 1};
    end
end

%% Display anatomical images
figure('name', 'Anatomic');
for i_eTE = 1:neTE
    subplot(1,neTE,i_eTE);
    imagesc(abs(Image_combined{i_eTE, 1}));caxis([0, 500]);
    title(sprintf('eTE = %.1f ms', eTE(i_eTE)));
    set(gca,'fontsize',15);
    colorbar;
    axis image;colormap gray;
    axis off;
end
saveas(gcf, [results_folder, filesep, 'Anatomic_all_eTE.fig']);

%% Display magnitude combined CD images
CD_avg_mag_combined = cell(neTE, 1);
for i_eTE = 1:neTE
    CD_avg_mag_combined{i_eTE} = zeros(Nrow, Ncol);
    
    for iSet = 1:nSet-1
        CD_avg_mag_combined{i_eTE} = CD_avg_mag_combined{i_eTE} + abs(CD_complex{i_eTE, iSet}).^2;
    end
    
    CD_avg_mag_combined{i_eTE} = sqrt(CD_avg_mag_combined{i_eTE});
end

for i_eTE = 1:neTE
    figure('name', sprintf('CD_avg_combined eTE = %.1f ms.fig', eTE(i_eTE)));
    imagesc(CD_avg_mag_combined{i_eTE});caxis(MagAxisSame);
    axis image;colormap gray;title(sprintf('eTE = %.1f ms', eTE(i_eTE)));
    colorbar;
    set(gca,'fontsize',15);
    axis off;
    saveas(gcf, [results_folder, filesep, sprintf('CD_avg_combined_eTE%.0f.fig', eTE(i_eTE))]);
end

save([results_folder, filesep, savefilename], 'CD_correct', 'CD', 'mag1', 'mag2', 'dphi', 'dphi_correct',...
    'CD_complex', 'im_orientation_PE', 'hct', 'eTE', 'T2fit_method', 'BW_static', 'BW_brain', 'rep2exclude');

CD4ROIdrawing = zeros(size(CD_complex{1, 1}));
for iSet = 1:nSet-1
    CD4ROIdrawing = CD4ROIdrawing + abs(CD_complex{1, iSet}).^2;
end
CD4ROIdrawing = sqrt(CD4ROIdrawing);

%% T2 and Yv map fitting

if bT2Yvmap && length(eTE(ind2fit)) > 1   % temp remove

    %% Draw vessel mask
    Im2disp = CD4ROIdrawing;

    fig_vessel_ROI = figure('name', 'Vessel mask');
    imshow(Im2disp, MagAxisSame);hold on;
    addToolbarExplorationButtons(fig_vessel_ROI); % this is important, otherwise using zoom-in during roipoly() is super annoying!
    
    if bLoadPreVesselMask
        [BW_vessel_filename, BW_vessel_pathname] = uigetfile('*.mat', 'Select preliminary vessel mask');
        load([BW_vessel_pathname, BW_vessel_filename], 'BW_vessel', 'x_vessel', 'y_vessel');
        
        NumVesselROI = length(BW_vessel);
        for iROI = 1:NumVesselROI
            plot(x_vessel{iROI}, y_vessel{iROI}, 'r', 'linewidth', 2);
        end
    else
%         NumVesselROI = str2double(inputdlg('Number of vessel ROI'));     
        NumVesselROI = 1;  % tem change
        BW_vessel = cell(NumVesselROI, 1);
        x_vessel = cell(NumVesselROI, 1);
        y_vessel = cell(NumVesselROI, 1);
        for iROI = 1:NumVesselROI
            title(sprintf('Draw vessel mask: #%d', iROI));set(gca,'fontsize',15);
            [BW_vessel{iROI}, x_vessel{iROI}, y_vessel{iROI}] = roipoly();
            plot(x_vessel{iROI}, y_vessel{iROI}, 'r', 'linewidth', 2);
        end
    end
    hold off;
    drawnow;
    zoom out;
    save([results_folder, filesep,'pre_vessel_mask.mat'], 'BW_vessel', 'x_vessel', 'y_vessel');
    saveas(gcf, [results_folder, filesep, 'pre_vessel_mask.fig']);

%%  replace above: use processed vessel mask
%     cd(pathname);
%     load('./aTRUPC_mid_sagittal_real_imag_Auto_Manual/pre_vessel_mask.mat'); % add
%     Im2disp = CD4ROIdrawing;
%     NumVesselROI = length(BW_vessel);
%%    
    mask_vessel_combined = zeros(size(BW_vessel{1}));
    for iROI = 1:NumVesselROI
        mask_vessel_combined = mask_vessel_combined | BW_vessel{iROI};
    end
    
    % Thresholding
    vessel_thresh = median(Im2disp(mask_vessel_combined))*vessel_thresh_ratio;
    mask_vessel_final = Im2disp.*mask_vessel_combined > vessel_thresh;
    figure('name','vessel mask'),imshow(mask_vessel_final);
    drawnow;
    saveas(gcf, [results_folder, filesep, 'final_vessel_mask.fig']);
    % save vessel mask
    save([results_folder, filesep, 'vesselmaskloc.mat'], 'mask_vessel_final');

    %% T2 and Yv map fitting in a sliding-window manner 
    T2map = zeros(Nrow, Ncol);
    Yvmap = zeros(Nrow, Ncol);
    DeltaR2map = zeros(Nrow, Ncol);
    
    if bT2Correction
        T2map_corrected = zeros(Nrow, Ncol);
        Yvmap_corrected = zeros(Nrow, Ncol);
    end
    
    % use all vessel pixels in the sliding window to fit for T2
    for iRow = 1:Nrow
        for iCol = 1:Ncol
            if mask_vessel_final(iRow, iCol) > 0
                RowMin = max(iRow-floor(T2map_window_size(1)/2),1);
                RowMax = min(iRow-floor(T2map_window_size(1)/2)+T2map_window_size(1)-1, Nrow);
                ColMin = max(iCol-floor(T2map_window_size(2)/2),1);
                ColMax = min(iCol-floor(T2map_window_size(2)/2)+T2map_window_size(2)-1, Ncol);
                mask_block = mask_vessel_final(RowMin:RowMax, ColMin:ColMax);
                NumPix = sum(mask_block(:));
                
                eTE_All = eTE(ind2fit);
                eTE_All = eTE_All(:);
                if strcmpi(T2fit_method, 'multiS0')
                    eTE_All = repmat(eTE_All, [NumPix*2*(nSet-1), 1]);
                end
                
                sig_window_real = zeros(length(ind2fit), NumPix*(nSet-1));
                sig_window_imag = zeros(size(sig_window_real));
                PixInd_start = 1;
                for iSet = 1:nSet-1
                    for iIM = 1:length(ind2fit)
                        image_block = CD_complex{ind2fit(iIM), iSet}(RowMin:RowMax, ColMin:ColMax);
                        sig_window_real(iIM, PixInd_start:PixInd_start+NumPix-1) = real(image_block(mask_block));
                        sig_window_imag(iIM, PixInd_start:PixInd_start+NumPix-1) = imag(image_block(mask_block));
                    end
                    PixInd_start = PixInd_start + NumPix;
                end
                sig_window = [sig_window_real, sig_window_imag];
                if strcmpi(T2fit_method, 'multiS0')
                    sig2fit = sig_window;
                    sig2fit = sig2fit(:);
                    init_S0 = sig_window(1,:);
                    init_beta = [init_S0(:); 1/50];
                    [beta_bar, resid_bar, jacob_bar] = nlinfit_hlu([eTE_All;length(eTE(ind2fit))], sig2fit, 'multi_monexp_model_nte', init_beta);
                elseif strcmpi(T2fit_method, 'avg')
                    [beta_bar, resid_bar, jacob_bar] = nlinfit_hlu(eTE_All, mean(sig_window, 2), 'MonoExp_model',[100,1/50]);
                end

                T2map(iRow, iCol) = 1./beta_bar(end);
                Yvmap(iRow, iCol) = 100*T2toY(T2map(iRow, iCol), hct);
                conintval1 = nlparci(beta_bar, resid_bar, jacob_bar); %95% confidence interval for estimates
                DeltaR2map(iRow, iCol) = 1000*(conintval1(end,2)-conintval1(end,1));
                
                if bT2Correction
                    T2map_corrected(iRow, iCol) = Fun_T2Cor(T2map(iRow, iCol));
                    Yvmap_corrected(iRow, iCol) = 100*T2toY(T2map_corrected(iRow, iCol), hct);
                end
            end
        end
    end
    
    %% Display T2 and Yv map
    T2caxis = [40,120];
    T2cmap = jet;
    ind_0 = ceil(size(T2cmap,1)*(0-T2caxis(1))/(T2caxis(2)-T2caxis(1))); % row index in the colormap for value 0
    ind_0 = min(max(ind_0, 1),size(T2cmap,1));
    T2cmap(ind_0,:) = [0,0,0]; % set color for value 0 to black
        
    Yvcaxis = [50,100];
    Yvcmap = jet;
    ind_0 = ceil(size(Yvcmap,1)*(0-Yvcaxis(1))/(Yvcaxis(2)-Yvcaxis(1))); % row index in the colormap for value 0
    ind_0 = min(max(ind_0, 1),size(Yvcmap,1));
    Yvcmap(ind_0,:) = [0,0,0]; % set color for value 0 to black
    
    if bT2Correction
        figure('name','T2 map corrected');
        imagesc(T2map_corrected);axis image;
        colorbar;set(gca,'fontsize',20,'fontweight','bold');
        caxis(T2caxis);
        colormap(T2cmap);
        axis off;
        title('T2 map corrected (ms)');
        saveas(gcf, [results_folder, filesep, 'T2map_corrected.fig']);

        figure('name','Yv map corrected');
        imagesc(Yvmap_corrected);axis image;
        colorbar;set(gca,'fontsize',20,'fontweight','bold');
        caxis(Yvcaxis);
        colormap(Yvcmap);
        axis off;
        title('Yv map corrected (%)');
        saveas(gcf, [results_folder, filesep, 'Yvmap_corrected.fig']);
    else
        figure('name','T2 map');
        imagesc(T2map);axis image;
        colorbar;set(gca,'fontsize',20,'fontweight','bold');
        caxis(T2caxis);
        colormap(T2cmap);
        axis off;
        title('T2 map (ms)');
        saveas(gcf, [results_folder, filesep, 'T2map.fig']);

        figure('name','Yv map');
        imagesc(Yvmap);axis image;
        colorbar;set(gca,'fontsize',20,'fontweight','bold');
        caxis(Yvcaxis);
        colormap(Yvcmap);
        axis off;
        title('Yv map (%)');
        saveas(gcf, [results_folder, filesep, 'Yvmap.fig']);
    end
    
    figure('name','Delta-R2 map');
    DeltaR2caxis = [0,5];
    DeltaR2cmap = hot;
    ind_0 = ceil(size(DeltaR2cmap,1)*(0-DeltaR2caxis(1))/(DeltaR2caxis(2)-DeltaR2caxis(1))); % row index in the colormap for value 0
    ind_0 = min(max(ind_0, 1),size(DeltaR2cmap,1));
    DeltaR2cmap(ind_0,:) = [0,0,0]; % set color for value 0 to black
    imagesc(DeltaR2map);axis image;
    colorbar;set(gca,'fontsize',20,'fontweight','bold');
    caxis(DeltaR2caxis);
    colormap(DeltaR2cmap);
    axis off;
    title('Delta-R2 (Hz)');
    saveas(gcf, [results_folder, filesep, 'DeltaR2map.fig']);
    
    save([results_folder, filesep, savefilename], 'T2map', 'Yvmap', 'DeltaR2map',...
        'BW_vessel', 'x_vessel', 'y_vessel', '-append');
    
    if bT2Correction
        save([results_folder, filesep, savefilename], 'T2map_corrected', 'Yvmap_corrected', '-append');
    end
end   % temp remove

%% Select the SSS, Straight Sinus, Deep Brain veins ROI manually
fig_ROI3 = figure('name','select three ROIs');
addToolbarExplorationButtons(fig_ROI3); % this is important, otherwise using zoom-in during roipoly() is super annoying!
imshow(mask_vessel_final);

BW = cell(1,1);
ROIx = cell(1,1);
ROIy = cell(1,1);
hold on;
for iROI = 1:3
    if iROI == 1
        title(sprintf('draw ROI#%d / %d: SSS', iROI, 3));
    elseif iROI == 2
        title(sprintf('draw ROI#%d / %d: Straight Sinus, include turning point', iROI, 3));
    else 
        title(sprintf('draw ROI#%d / %d: Deep Brain veins', iROI, 3));    
    end    
    [BW{iROI}, ROIx{iROI}, ROIy{iROI}] = roipoly();
    plot(ROIx{iROI},ROIy{iROI},'r','linewidth',2);
    text(mean(ROIx{iROI}), mean(ROIy{iROI}), sprintf('%d', iROI), 'color', [0.929, 0.694, 0.125],...
        'fontsize', 20, 'fontweight', 'bold');
end
hold off;

title(sprintf('ROI of three regions'));set(gca,'fontsize',15);
zoom out;

%% 20230627  bwskel: SSS

% Based on bwskel & region-growing alike method, three ROI selection

BW_SSS = BW{1}.* mask_vessel_final;
figure,imshow(BW_SSS);title('SSS mask');
branchlimit = 20;
[skelb_SSS,skel_SSS, coor_AP, coor_disp] = sortedskel(BW_SSS,1,branchlimit);

% suggest: Adult: physical length: 30mm
% pixel num = 30/(200/256) = 38.4

l_ROI = Physcdist_1/voxelsize(1); % pixel length of ROI

coor_cumudis = cumsum(coor_disp);
[~,ind] = min(abs(coor_cumudis - l_ROI));
coor_A(1,:)=coor_AP(1,:);coor_P(1,:)=coor_AP(ind,:); % ROI 1

[~,ind2] = min(abs(coor_cumudis - coor_cumudis(end,:)/2));
dist = (coor_cumudis(ind2)-l_ROI/2);
[~,ind3] = min(abs(coor_cumudis - dist));
coor_A(2,:)=coor_AP(ind3,:);
coor2_disp=coor_disp;coor2_disp(1:ind2,:)=[];
coor2_cumudis = cumsum(coor2_disp);
[~,ind4] = min(abs(coor2_cumudis - l_ROI/2));
coor_P(2,:)=coor_AP(ind4+ind2,:);  % ROI 2

coor3_disp=flip(coor_disp);
coor3_cumudis = cumsum(coor3_disp);
[~,ind] = min(abs(coor3_cumudis - l_ROI));
coor_A(3,:)=coor_AP(length(coor_disp)-ind,:);coor_P(3,:)=coor_AP(end,:); % ROI 3

BW_ROI = {};
figure, imagesc(BW_SSS);axis off;axis image;colormap gray;hold on;
plot(coor_AP(:,2),coor_AP(:,1));
for ROI_i = 1:3
    radii = 0.5*norm(coor_A(ROI_i,:)-coor_P(ROI_i,:));  % ROI 1 radii
    ROI_center = round(0.5*(coor_A(ROI_i,:)+coor_P(ROI_i,:)));
    circle_mask = createCirclesMask(size(BW_SSS),flip(ROI_center),radii);
    BW_ROI{ROI_i,1}=BW_SSS.*circle_mask;
%     figure,imshow([BW_SSS,BW_ROI{ROI_i,1}]);
    plot([coor_A(ROI_i,2);coor_P(ROI_i,2)],[coor_A(ROI_i,1);coor_P(ROI_i,1)],'b*');
end
hold off;

%% 20230707 bwskel: SS:

BW_SS = BW{2}.* mask_vessel_final; figure,imshow(BW_SS);title('SS mask');

branchlimit = 5;
[skelb_SS,skel_SS, coor_AP, coor_disp] = sortedskel(BW_SS,1,branchlimit);

coor_cumudis = cumsum(coor_disp);

[~,ind] = min(abs(coor_cumudis - coor_cumudis(end,:)/2));
% ind = round(size(coor_AP,1)/2);

dist = (coor_cumudis(ind)-l_ROI/2);
[~,ind2] = min(abs(coor_cumudis - dist));
coor_A(4,:)=coor_AP(ind2,:);

coor2_disp=coor_disp;coor2_disp(1:ind,:)=[];
coor2_cumudis = cumsum(coor2_disp);
[~,ind3] = min(abs(coor2_cumudis - l_ROI/2));
coor_P(4,:)=coor_AP(ind3+ind,:); % ROI 2

radii = 0.5*norm(coor_A(4,:)-coor_P(4,:));  % ROI 1 radii
ROI_center = round(0.5*(coor_A(4,:)+coor_P(4,:)));
circle_mask = createCirclesMask(size(BW_SS),flip(ROI_center),radii);
BW_ROI{4,1}=BW_SS.*circle_mask;
% figure,imshow([BW_SS,BW_ROI{4,1}]);

%verify
figure, imagesc(BW_SS);axis off;axis image;colormap gray;hold on;
plot(coor_AP(:,2),coor_AP(:,1),[coor_A(4,2);coor_P(4,2)],[coor_A(4,1);coor_P(4,1)],'b*')
plot(coor_AP(ind,2),coor_AP(ind,1),'r*');hold off;

%% 20230711 bwskel: GV/ICV:

BW_DB0 = BW{3}.* mask_vessel_final;
% figure,imshow(BW_DB0);title('DB mask');
BW_DB = logical(imfill(BW_DB0,'holes'));
out = bwskel(BW_DB);
out2 = bwskel(BW_DB,'MinBranchLength',5); % the size of the disk needs modification

% Original ske
figure, imagesc(BW_DB0);axis off;axis image;colormap gray;hold on;
[row,col] = find(out2==1);
plot(col,row,'r.');
hold off;

% 20230714 update: curve fitting on ske to extend it to object border
p = polyfit(col,row,5);
xcol0 = linspace(min(col)-2, max(col)+2)'; % extend two ends
yrow0 = polyval(p,xcol0);
figure, imagesc(BW_DB0);axis off;axis image;colormap gray;hold on; % ori and extended ske
plot(col,row,'b.','LineWidth',2);
plot(xcol0,yrow0,'r.','LineWidth',2);
hold off;
% find unique pair
ind = unique(sub2ind(size(BW_DB0),round(yrow0),round(xcol0)));
skefit = zeros(size(BW_DB0));
skefit(ind)=1; skefit = skefit.*BW_DB;
figure, imagesc(BW_DB0);axis off;axis image;colormap gray;hold on;
[yrow,xcol] = find(skefit); %pixelized curve fitting (--extended skeleton)
plot(xcol,yrow,'r.');
hold off;
[~, ~,coor_AP,coor_disp] = sortedskel(skefit, 0);%--extended and sorted skeleton

% find anterior & posterior point on DB
coor_A(6,:) = coor_AP(1,:);
coor_P(5,:) = coor_AP(end,:);

% find extreme point: local minima(treat it as the turning point)
figure, imagesc(BW_DB0);axis off;axis image;colormap gray;hold on;
TF = islocalmin(size(BW_DB0,1)-yrow0);
plot(xcol0,yrow0,xcol0(TF),yrow0(TF),'r*')
hold off
title(sprintf('Local minima'));set(gca,'fontsize',15);
% saveas(gcf, [results_folder, filesep, 'ROI.fig']);

% find the point seperating GV & ICV (rightmost local minima)
cnpoint = [yrow0(TF),xcol0(TF)];
ppoint = cnpoint(end,:);

% find GV end
dist = coor_AP-ceil(ppoint);
[~,ind2] = min(dist(:,1).^2+dist(:,2).^2);
coor_A(5,:) =coor_AP(ind2,:);
coor_cumudis = cumsum(coor_disp);

% from turning point, move 10mm (1cm) towards the end of ICV along the skeleton
% 12.8 *200/256 =   10mm
l_interval = Physcdist_2/voxelsize(1); % pixel length of two ROI interval

dist = coor_cumudis(ind2) - l_interval;
[~,ind3] = min(abs(coor_cumudis - dist));
coor_P(6,:)=coor_AP(ind3,:); % find ICV end

% find radii
ICVr = norm(coor_A(6,:) - coor_P(6,:));
GVr = norm(coor_A(5,:) - coor_P(5,:));
 
radii = ceil(GVr);
% ROI5 center of circle
ROI_center = coor_P(5,:);
circle_mask = createCirclesMask(size(BW_DB0),flip(ROI_center),radii);
BW_ROI{5,1}=BW_DB0.*circle_mask;

radii = ceil(ICVr);
% ROI6 center of circle
ROI_center = coor_A(6,:);
circle_mask = createCirclesMask(size(BW_DB0),flip(ROI_center),radii);
BW_ROI{6,1}=BW_DB0.*circle_mask;

% %verify
% figure, imagesc(BW_DB0);axis off;axis image;colormap gray;hold on;
% plot(xcol0,yrow0,'LineWidth',1.5);hold on
% plot(xcol0(TF),yrow0(TF),'r*')
% hold on
% plot([coor_A(6,2);coor_P(6,2)],[coor_A(6,1);coor_P(6,1)],'b*')


% plot 6 ROIs
figure, imagesc(mask_vessel_final);axis off;axis image;colormap gray;hold on;
for i=1:length(BW_ROI)
    [ROIy{i},ROIx{i}] = find(BW_ROI{i}==1);    
    plot(ROIx{i},ROIy{i},'.b','linewidth',2);
end
hold off;
title(sprintf('ROI for T2 and Yv'));set(gca,'fontsize',15);
save([results_folder, filesep, 'ROI2fit.mat'], 'BW_ROI','ROIx','ROIy');
saveas(gcf, [results_folder, filesep, 'ROI.fig']);


%% Select the pixels with highest magnitude of complex difference of either flow-encoding -- for six ROIs
ROIBW_MaxPixs = cell(NumROI, 1);
nte = length(ind2fit);

CD2pickpixel = zeros(Nrow, Ncol, nSet-1);
for iSet = 1:nSet-1
    CD2pickpixel(:,:,iSet) = abs(CD_complex{1, iSet});
end

for iROI = 1:NumROI
    BW_all = reshape(BW_ROI{iROI}, [Nrow, Ncol, 1]);
    BW_all = repmat(BW_all, [1,1,nSet-1]);
    if bUseMaxPixsInROI
        ROIBW_MaxPixs{iROI} = zeros(size(CD2pickpixel));
        cntIm = CD2pickpixel.*BW_all;
        [~,ind] = sort(cntIm(:), 'descend');
        % 20230710 update: in case the pixels are fewer than the threshold
        if NumPixInROI(iROI) > length(find(BW_all==1))
            ROIBW_MaxPixs{iROI} = logical(BW_all);
            NumPixInROI(iROI) = sum(BW_all(:));
        else
        ind4avg = ind(1:NumPixInROI(iROI));
        ROIBW_MaxPixs{iROI}(ind4avg) = 1;
        ROIBW_MaxPixs{iROI} = logical(ROIBW_MaxPixs{iROI});
        end
    else
        ROIBW_MaxPixs{iROI} = logical(BW_all);
        NumPixInROI(iROI) = sum(BW_all(:));
    end
end

% Display the final mask
figure('name', 'Final mask ROI');
for iSet = 1:nSet-1
    subplot(1, nSet-1, iSet);
    imagesc(abs(CD_complex{1, iSet}), MagAxisSame);axis off;axis image;colormap gray;
    hold on;
    for iROI = 1:NumROI
        plot(ROIx{iROI},ROIy{iROI},'.b');
        [DotY, DotX] = find(ROIBW_MaxPixs{iROI}(:,:,iSet) == 1);
        plot(DotX, DotY, '.r');
    end
    hold off;
end
saveas(gcf, [results_folder, filesep, 'Final_mask_ROI.fig']);

%% Do T2 and Yv fitting for real and imagninary parts
% NumROI = 3;
sig_ROI = cell(NumROI, 1);
sig_ROI_mag = cell(NumROI, 1);
sig_ROI_avg = zeros(length(ind2fit), NumROI, 2);
sig_ROI_mag_avg = zeros(length(ind2fit), NumROI);
T2_roi_AllRep = zeros(1, NumROI);
Yv_roi_AllRep = zeros(size(T2_roi_AllRep));
deltaR2_roi_AllRep = zeros(size(T2_roi_AllRep));

for iROI = 1:NumROI
    eTE_All = eTE(ind2fit);
    eTE_All = eTE_All(:);
    if strcmpi(T2fit_method, 'multiS0')
        eTE_All = repmat(eTE_All, [NumPixInROI(iROI)*2, 1]);
    end

    PixInd_start = 1;
    sig_ROI_real = zeros(length(ind2fit), NumPixInROI(iROI));
    sig_ROI_imag = zeros(size(sig_ROI_real));
    for iSet = 1:nSet-1
        NumPixInSet = sum(sum(ROIBW_MaxPixs{iROI}(:,:,iSet)));
        for iIM = 1:length(ind2fit)
            cntIm_real = real(CD_complex{ind2fit(iIM), iSet});
            cntIm_imag = imag(CD_complex{ind2fit(iIM), iSet});
            sig_ROI_real(iIM, PixInd_start:PixInd_start+NumPixInSet-1) = cntIm_real(ROIBW_MaxPixs{iROI}(:,:,iSet));
            sig_ROI_imag(iIM, PixInd_start:PixInd_start+NumPixInSet-1) = cntIm_imag(ROIBW_MaxPixs{iROI}(:,:,iSet));
        end
        PixInd_start = PixInd_start + NumPixInSet;
    end
    sig_ROI{iROI} = [sig_ROI_real, sig_ROI_imag];
	sig_ROI_mag{iROI} = abs(sig_ROI_real + sqrt(-1)*sig_ROI_imag);
    
    for iFitting = 1:2 % 1=real, 2=imaginary
        for iIM = 1:length(ind2fit)
            sig_ROI_avg(iIM, iROI, iFitting) = mean(sig_ROI{iROI}(iIM, (iFitting - 1)*NumPixInROI(iROI)+1:iFitting*NumPixInROI(iROI)));
        end
    end
    
	for iIM = 1:length(ind2fit)
        sig_ROI_mag_avg(iIM, iROI) = mean(sig_ROI_mag{iROI}(iIM, :));
    end
    % ROI T2 and Yv fitting
    if length(eTE(ind2fit)) > 1
        if strcmpi(T2fit_method, 'multiS0')
            sig2fit = sig_ROI{iROI};
            sig2fit = sig2fit(:);
            init_S0 = sig_ROI{iROI}(1,:);
            init_beta = [init_S0(:); 1/50];
            [beta_roi, resid_roi, jacob_roi] = nlinfit_hlu([eTE_All;length(eTE(ind2fit))], sig2fit, 'multi_monexp_model_nte', init_beta);
        elseif strcmpi(T2fit_method, 'avg')
            [beta_roi, resid_roi, jacob_roi] = nlinfit_hlu(eTE_All, mean(sig_ROI{iROI}, 2), 'MonoExp_model',[100,1/50]);
        end

        T2_roi_AllRep(iROI) = 1./beta_roi(end);
        Yv_roi_AllRep(iROI) = 100*T2toY(T2_roi_AllRep(iROI), hct);
        conintval1 = nlparci(beta_roi, resid_roi, jacob_roi); %95% confidence interval for estimates
        deltaR2_roi_AllRep(iROI) = 1000*(conintval1(end,2)-conintval1(end,1));
    end
end

% T2 correction
if length(eTE(ind2fit)) > 1
    if bT2Correction
        T2_roi_AllRep_corrected = zeros(size(T2_roi_AllRep));
        Yv_roi_AllRep_corrected = zeros(size(T2_roi_AllRep));
        for iROI = 1:NumROI
            T2_roi_AllRep_corrected(iROI) = Fun_T2Cor(T2_roi_AllRep(iROI));
            Yv_roi_AllRep_corrected(iROI) = 100*T2toY(T2_roi_AllRep_corrected(iROI), hct);
        end
    end
end

save([results_folder, filesep, savefilename], 'sig_ROI_mag', 'sig_ROI_mag_avg',...
    'BW', 'ROIx', 'ROIy', 'sig_ROI', 'sig_ROI_avg', 'NumPixInROI', 'bUseMaxPixsInROI', 'ind2fit', '-append');

if length(eTE(ind2fit)) > 1
    save([results_folder, filesep, savefilename], 'T2_roi_AllRep', 'Yv_roi_AllRep', 'deltaR2_roi_AllRep', '-append');
    
    if bT2Correction
        save([results_folder, filesep, savefilename], 'T2_roi_AllRep_corrected', 'Yv_roi_AllRep_corrected', '-append');
    end
end

%% Write ROI results into .csv file
if length(eTE(ind2fit)) > 1
    if bT2Correction
        rowTitle = {'ROI'; 'T2 corrected (ms)'; 'Delta-R2 (1/s)'; 'Yv corrected (%)'};
        output_mat = [T2_roi_AllRep_corrected; deltaR2_roi_AllRep; Yv_roi_AllRep_corrected];
    else
        rowTitle = {'ROI'; 'T2 (ms)'; 'Delta-R2 (1/s)'; 'Yv (%)'};
        output_mat = [T2_roi_AllRep; deltaR2_roi_AllRep; Yv_roi_AllRep];
    end
end
output_mat = [[1:NumROI]; output_mat];


fresult=[results_folder, filesep, 'ROI_results.csv'];
varNames = cell(NumROI, 1);
for iROI = 1:NumROI
    varNames{iROI} = sprintf('ROI%d', iROI);
end
T = table(output_mat, 'RowNames', rowTitle);
writetable(T, fresult, 'WriteRowNames', true, 'WriteVariableNames', false);

toc