function [CD_correct, dphi_correct, phi_fit, cc_fit, mask_static_refined] = CDPhaCor_Iter_2D_V3(mag1, mag2, dphi, Z1Z2product, mask_static, mask_brain, deg, Niter, bBatch)
% MODIFY ON 20230807: considering the batch condition

% this function correct the phase errors in complex difference iteratively,
% based on user input initial static tissue mask, and do the correction 
% either separately or for multiple images at a time;
% the static tissue mask is refined during the iterations, based on the
% results of the last iteration;
% the spatially dependent phase errors are fitted separately for each image
% or using data from all images conjunctively

% Input:
% mag1, mag2: cell arrays, combined magnitude images with bipolar gradient
% dphi: cell array, combined phase difference
% Z1Z2product: cell array, combined Z2*conj(Z1)
% mask_static: 2D, initial mask for static tissue
% mask_brain: brain mask
% deg: scalar, degree of polynomial to fit the phase of static tissue
% Niter: number of iterations
% bBatch: if true, the spatially dependent phase errors are fitted using
% data from all images conjunctively; if false, phase errors are fitted
% separately for each image
% ind:when batch, which rep to use 

% Output:
% CD_correct: cell array, corrected CD
% dphi_correct: cell array, corrected phase
% phi_fit: 2D matrix, fitted phase of static tissues
% cc_fit: 1D array, fitting coefficients
% mask_static_refined: final static tissue mask

% Dengrong Jiang, JHU BME, June 2023

if nargin < 9
    bBatch = 0;
end
nRep = length(mag1);
nData = numel(mag1);
mask_static_refined = mask_static;

% phase correction
CD_correct = cell(size(mag1));
dphi_correct = cell(size(mag1));

if ~bBatch % phase errors are fitted separately for each image
    cc_fit = cell(size(mag1));
    phi_fit = cell(size(mag1));
end

imsize = size(dphi{1});
if deg == 2 % quadratic fitting
    if numel(imsize) == 2 % 2D
        [dc, x, y, xy, x2, y2] = Calcfield(imsize);
        for iIter = 1:Niter
            if bBatch % phase phase errors are fitted using data from all images conjunctively;
                Atotal = [];
                rhs = []; % right hand side of the normal equation
                % construct system matrix using all images conjuctively
                for iData = 1:nData
                    index = find(mask_static_refined>0);
                    A1 = dc(index);
                    A2 = x(index);
                    A3 = y(index);
                    A4 = xy(index);
                    A5 = x2(index);
                    A6 = y2(index);
                    A=[A1 A2 A3 A4 A5 A6];
                    Atotal = [Atotal; A];
                    rhs = [rhs; dphi{iData}(index)];
                end
                % solve the normal equation by pseudo-inverse
                cc_fit=(Atotal.'*Atotal)\(Atotal.')*rhs; % coefficients
                phi_fit = cc_fit(1)*dc + cc_fit(2)*x + cc_fit(3)*y + cc_fit(4)*xy + cc_fit(5)*x2 + cc_fit(6)*y2;

                % apply phase correction to all images            
                for iData = 1:nData
                    dphi_correct{iData} = argument(Z1Z2product{iData}.*exp(-1*sqrt(-1)*phi_fit));
                    % In theory, a^2+b^2-2*a*b*cos(theta) >= 0, however, if
                    % you use single precision, in some cases, this can be
                    % negative, making CD imaginary. To improve robustness,
                    % use a zero cut-off
                    CD_correct{iData} = sqrt(max(0, (mag1{iData}.^2+mag2{iData}.^2-2*mag1{iData}.*mag2{iData}.*cos(dphi_correct{iData}))));
                end


            else % phase errors are fitted separately for each image
                for iData = 1:nData
                    % construct system matrix
                    index = find(mask_static_refined>0);
                    A1 = dc(index);
                    A2 = x(index);
                    A3 = y(index);
                    A4 = xy(index);
                    A5 = x2(index);
                    A6 = y2(index);
                    Atotal=[A1 A2 A3 A4 A5 A6];
                    rhs = dphi{iData}(index);
                    % solve the normal equation by pseudo-inverse
                    cc_fit{iData}=(Atotal.'*Atotal)\(Atotal.')*rhs; % coefficients
                    phi_fit{iData} = cc_fit{iData}(1)*dc + cc_fit{iData}(2)*x + cc_fit{iData}(3)*y + cc_fit{iData}(4)*xy + cc_fit{iData}(5)*x2 + cc_fit{iData}(6)*y2;
                    % apply phase correction
                    dphi_correct{iData} = argument(Z1Z2product{iData}.*exp(-1*sqrt(-1)*phi_fit{iData}));
                    % In theory, a^2+b^2-2*a*b*cos(theta) >= 0, however, if
                    % you use single precision, in some cases, this can be
                    % negative, making CD imaginary. To improve robustness,
                    % use a zero cut-off
                    CD_correct{iData} = sqrt(max(0, (mag1{iData}.^2+mag2{iData}.^2-2*mag1{iData}.*mag2{iData}.*cos(dphi_correct{iData}))));
                end
            end

            if iIter ==1
                % Use the dynamic that has the least L1 norm within the brain:
                for iRep = 1:nRep
                    CD_corabs{iRep,1} =  sqrt(CD_correct{iRep,1}.^2);
                    L1normCD(iRep,1) = sum(CD_corabs{iRep,1}.*mask_brain,'all');
                end
                [~,ind_L1min] = min(L1normCD);
            end

            % Use the dynamic that has the least L1 norm within the brain:
            
           
            % refine the static tissue mask
            if iIter < Niter
                CD2mask = CD_correct{ind_L1min,1} .* mask_brain;
%                 CD2mask = CD_correct{1} .* mask_brain;
                CD2mask_norm = CD2mask./max(CD2mask(:));
                vessel_thresh = graythresh(CD2mask_norm(logical(mask_brain)));
                vessel_mask = CD2mask_norm > vessel_thresh;

                figure, imshow(vessel_mask, []);

                vessel_mask = bwareaopen(vessel_mask, 2); % remove single dots
                dilate = 10;
                se = strel('cube', dilate);
                vessel_mask_dilate = imdilate(vessel_mask,se) .* mask_brain;
                mask_static_refined = mask_brain - vessel_mask_dilate;

            end
        end
    end
end




