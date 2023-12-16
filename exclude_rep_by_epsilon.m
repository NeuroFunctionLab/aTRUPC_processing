function epsilon = exclude_rep_by_epsilon(results_folder)
% 20230808 Calculate epsilon for repetitions

load([results_folder, filesep, 'varforrep.mat']);
%% Realign across repetitions
% repetition exclusion based on epsilon
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
else 
    % realign to the 1 st repetion
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

%% Adaptive combine multi-coil complex images
Image_combined = cell(nRep, neTE, nSet);
weights = cell(nRep,1);
% use the eTE0 flow-compensated images to estimate the weights
for iRep = 1:nRep
    [Image_combined{iRep,1,1}, weights{iRep,1}, ~] = CoilAdaptiveCombine_V2(Image_coil{iRep, 1, 1}, 1, 9, 1, eye(nCh));
end
% use same weights for each repetition
for iRep = 1:nRep
    for i_eTE = 1:neTE
        for iSet = 1:nSet
            if ~(i_eTE == 1 && iSet == 1)
                Image_combined{iRep, i_eTE, iSet} = sum(Image_coil{iRep, i_eTE, iSet} .* conj(weights{iRep,1}), 3);
            end
        end
    end
end

%% Segmentation of static tissue mask
CD_avg_mag_combined = cell(nRep, neTE, 1);

% Calculate complex difference for subsequent eddy current correction
CD = cell(nRep, neTE, nSet-1);
mag1 = cell(nRep, neTE, 1);
mag2 = cell(nRep, neTE, nSet-1);
dphi = cell(nRep, neTE, nSet-1);
Z1Z2product = cell(nRep, neTE, nSet-1);
L1normCD = zeros(nRep,1);
CD2mask= cell(nRep,1);

for iRep = 1:nRep
    for i_eTE = 1:neTE
        mag1{iRep,i_eTE} = abs(Image_combined{iRep, i_eTE, 1});
        for iSet = 1:nSet-1
            mag2{iRep, i_eTE, iSet} = abs(Image_combined{iRep, i_eTE, iSet+1});
            Z1Z2product{iRep, i_eTE, iSet} = Image_combined{iRep, i_eTE, iSet+1}.*conj(Image_combined{iRep, i_eTE, 1});
            dphi{iRep, i_eTE, iSet} = argument(Z1Z2product{iRep, i_eTE, iSet});
            CD{iRep, i_eTE, iSet} = abs(Image_combined{iRep, i_eTE, iSet+1} - Image_combined{iRep, i_eTE, 1});
        end
    end
end 

%% Pre-segmentation
% Use anatomical image:
AI = abs(Image_combined{1, 1, 1});
AI_2 = AI./max(AI(:)); 

k1=AI_2 > graythresh(AI_2);
SE=strel('disk',7,4);
k2=imopen(k1,SE); 
b=bwlabel(k2);% apply connected component analysis.
[m,n]=size(b);
index=b(m.*0.5,n.*0.5);
b(b~=index)=0;
BW_brain=imclose(b,strel('disk',18)); 
C2 = labeloverlay(AI_2,BW_brain);

% pre-tissue selection: use the dynamic that has the least L1 norm of the whole image 
for iRep = 1:nRep
    CD2mask{iRep,1} =  sqrt(CD{iRep,1,1}.^2);
    L1normCD(iRep,1) = sum(CD2mask{iRep,1}.*BW_brain,'all');
end
[~,ind_L1min] = min(L1normCD);
figure('name', 'pre-select the frame with least L1 norm');   
for iRep = 1:nRep
    subplot(1, nRep, iRep);
    imagesc(CD2mask{iRep, 1});caxis(MagAxisSame);
    axis image;colormap gray;
    set(gca,'fontsize',15);axis off;
    if iRep == ind_L1min
        title(sprintf('Rep #%d, L1 norm = %.3e', iRep, L1normCD(iRep,1)),'Color', 'b');
    else
        title(sprintf('Rep #%d, L1 norm = %.3e', iRep, L1normCD(iRep,1)));
    end    
end 

CD2mask_min = CD2mask{ind_L1min, 1};
% threshold to segment out the vessel
CD2mask_brain = CD2mask_min.*BW_brain; 
CD2mask_brain_2 = CD2mask_brain./max(CD2mask_brain(:)); 
BW_vessel = CD2mask_brain_2 > graythresh(CD2mask_brain_2);
CD_image = CD2mask_min./max(CD2mask_min(:));
C1 = labeloverlay(CD_image,BW_vessel);
% substract to get static tissue image.
BW_static = BW_brain - BW_vessel;
C = labeloverlay(CD_image,BW_static);

%% Automatic eddy current correction
CD_correct = cell(nRep, neTE, nSet-1);
dphi_correct = cell(nRep, neTE, nSet-1);
Niter = 5;
bBatch = 1;
for iSet = 1:nSet-1
    [CD_correct, dphi_correct, phi_fit, cc_fit, mask_static_refined] = CDPhaCor_Iter_2D_V3...
        (mag1, mag2(:, :, iSet), dphi(:, :, iSet), Z1Z2product(:, :, iSet), BW_static, BW_brain, deg_poly_eddycurrent, Niter, bBatch);
end

C3 = labeloverlay(CD_image,mask_static_refined);
figure,imshow(C3);title('Static tissue ROI for calculating epsilon');
saveas(gcf, [results_folder, filesep, 'STforepsilon.fig']);

%%  Calculate epsilon 
epsilon = zeros(nRep,1);
for iRep = 1:nRep
    CD_complex = cell(neTE, nSet-1);
    for i_eTE = 1:neTE
        for iSet = 1:nSet-1
            CD_complex{i_eTE, iSet} = Image_combined{iRep, i_eTE, iSet+1}.*exp(-sqrt(-1).*phi_fit) - Image_combined{iRep, i_eTE, 1};
        end
    end
    
    for i_eTE = 1:neTE
        CD_avg_mag_combined{iRep, i_eTE} = zeros(Nrow, Ncol);
        
        for iSet = 1:nSet-1
            CD_avg_mag_combined{iRep, i_eTE} = CD_avg_mag_combined{iRep, i_eTE} + abs(CD_complex{i_eTE, iSet}).^2;
        end
        
        CD_avg_mag_combined{iRep, i_eTE} = sqrt(CD_avg_mag_combined{iRep, i_eTE});
    end
    epsilon(iRep) = mean(abs(mask_static_refined.*CD_avg_mag_combined{iRep,1}),'all') ./ mean(abs(mask_static_refined.*Image_combined{iRep,1,1}),'all')
end 
%% Display epsilon with CD images

figure('name', 'Inspect repetitions', 'WindowState', 'maximized');
for iRep = 1:nRep
    subplot(2, nRep, iRep);
    imagesc(CD_avg_mag_combined{iRep, 1});caxis(MagAxisSame);
    axis image;colormap gray;
    set(gca,'fontsize',15);axis off;

    title(sprintf('Rep #%d, E = %.4f', iRep, epsilon(iRep)));
%     title(' $\varepsilon$ ','interpreter','latex')
    subplot(2, nRep, iRep+nRep);
    imagesc(abs(Image_combined{iRep,1,1}));caxis([0, 500]);
    axis image;colormap gray;axis off;
end
saveas(gcf, [results_folder, filesep, 'inspect_repetitions_epsilon.fig']); % add this save file...

delete([results_folder, filesep, '*.hdr']);
delete([results_folder, filesep, '*.img']);
