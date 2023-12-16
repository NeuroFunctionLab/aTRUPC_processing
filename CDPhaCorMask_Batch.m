function [CD_correct, dphi_correct, phi_fit, cc_fit, mask_reg] = CDPhaCorMask_Batch(mag1, mag2, dphi, Z1Z2product, mask_static, deg)
% this function correct the phase errors in complex difference, based on
% user input static tissue mask, and do the correction for multiple images
% the spatially dependent phase errors are fitted using data from all 
% images conjunctively

% Input:
% mag1, mag2: cell arrays, combined magnitude images with bipolar gradient
% dphi: cell array, combined phase difference
% Z1Z2product: cell array, combined Z2*conj(Z1)
% mask_static: 2D matrix, mask for static tissue, assumed to be drawn on
% the complex difference image corresponding to the 1st image in mag1
% deg: scalar, degree of polynomial to fit the phase of static tissue

% Output:
% CD_correct: cell array, corrected CD
% dphi_correct: cell array, corrected phase
% phi_fit: 2D matrix, fitted phase of static tissues
% cc_fit: 1D array, fitting coefficients
% mask_reg: cell array, mask_static registered to each individual image

% Dengrong Jiang, JHU BME, Dec 2017
nData = length(mag1);

mask_reg = mask_static; % use the 1st image as the reference image

imsize = size(dphi{1});
if deg == 2 % quadratic fitting
    if numel(imsize) == 2 % 2D
        [dc, x, y, xy, x2, y2] = Calcfield(imsize);
        Atotal = [];
        rhs = []; % right hand side of the normal equation
        % construct system matrix using all images conjuctively
        for iData = 1:nData
            index = find(mask_reg > 0);
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
        CD_correct = cell(size(mag1));
        dphi_correct = cell(size(mag1));
        for iData = 1:nData
            dphi_correct{iData} = argument(Z1Z2product{iData}.*exp(-1*sqrt(-1)*phi_fit));
            CD_correct{iData} = sqrt(mag1{iData}.^2+mag2{iData}.^2-2*mag1{iData}.*mag2{iData}.*cos(dphi_correct{iData}));
        end
    end
end




