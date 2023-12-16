function [Im, wtCoarse, wtFull] = CoilAdaptiveCombine_V2(Im_coil, blocksize, trainingsize, bInterp, Rn)
% 20230615 update: line 98 changed: delete sqrt

% Input:
% Im_coil: Nrow x Ncol x Ncoil 256 256 32
% blocksize: size of the image block
% trainingsize: size of the training region
% bInterp: if true, interpolate the weights to full image matrix, instead
% of using the same weights for a whole image block
% Rn: noise correlation matrix

% This function performs adaptive combine for multi-coil images, based on
% D. Walsh et al. MRM 43:682-690 (2000)
% Phase correction and interpolation is based on Griswold et al. ISMRM 2002 p2410
% Reference codes: https://github.com/andyschwarzl/gpuNUFFT/blob/master/matlab/demo/utils/adapt_array_2d.m
% Dengrong Jiang, JHU BME Dec. 2017

[Nrow, Ncol, Ncoil] = size(Im_coil);

% if Rn is not obtained externally, draw noise ROI and compute Rn within the ROI
if nargin < 5
    sos = sqrt(sum(Im_coil.*conj(Im_coil), 3));
    figure,imshow(sos, []);title('Please draw noise ROI', 'fontsize', 15);
    [NoiseMask] = roipoly();
    close;
    % construct noise correlation matrix
    % Direct implementation
    % [Row_noise, Col_noise] = find(NoiseMask == 1);
    % Rn = zeros(Ncoil, Ncoil); % noise correlation matrix
    % for iPx = 1:length(Row_noise)
    %     cntPx = Im_coil(Row_noise(iPx), Col_noise(iPx), :);
    %     cntPx = cntPx(:);
    %     Rn = Rn + cntPx * cntPx';
    % end

    % Vectorized implementation, equivalent to direct implementation
    Ind_noise = find(NoiseMask == 1);
    Rn_pre = zeros(length(Ind_noise), Ncoil);
    for iCoil = 1:Ncoil
        cntIm = Im_coil(:,:,iCoil);
        Rn_pre(:, iCoil) = cntIm(Ind_noise);
    end
    Rn = Rn_pre.'*conj(Rn_pre);
    Rn = Rn/norm(Rn, 'fro')*sqrt(Ncoil); % normalize Rn towards identity matrix
end

if nargin < 4
    bInterp = 0;
end

% find coil with max intensity for phase correction
% Based on Griswold et al. ISMRM 2002 p2410
[~, maxCoil] = max(sum(sum(abs(Im_coil),1),2));

% inverse of noise correlation matrix
Rn_inv = inv(Rn);

%% block-by-block computation
Im = zeros(Nrow, Ncol);
wtCoarse = zeros(length([1:blocksize:Nrow]), length([1:blocksize:Ncol]), Ncoil);

for iRow = blocksize:blocksize:Nrow
    for iCol = blocksize:blocksize:Ncol
        % [iRow, iCol] is the center of the image block and training block
        % compute signal correlation matrix
        % Direct implementation
%         TrainRowMin = max(iRow-floor(trainingsize/2), 1);
%         TrainRowMax = min(iRow-floor(trainingsize/2)+trainingsize-1, Nrow);
%         TrainColMin = max(iCol-floor(trainingsize/2), 1);
%         TrainColMax = min(iCol-floor(trainingsize/2)+trainingsize-1, Ncol);
%         Rs = zeros(Ncoil, Ncoil); % signal correlation matrix
%         for iTrainRow = TrainRowMin:TrainRowMax
%             for iTrainCol = TrainColMin:TrainColMax
%                 cntPx = Im_coil(iTrainRow, iTrainCol, :);
%                 cntPx = cntPx(:);
%                 Rs = Rs + cntPx * cntPx';
%             end
%         end

        % Vectorized implementation, equivalent to direct implementation
        TrainRowMin = max(iRow-floor(trainingsize/2), 1);
        TrainRowMax = min(iRow-floor(trainingsize/2)+trainingsize-1, Nrow);
        TrainColMin = max(iCol-floor(trainingsize/2), 1);
        TrainColMax = min(iCol-floor(trainingsize/2)+trainingsize-1, Ncol);
        Rs_pre = Im_coil(TrainRowMin:TrainRowMax, TrainColMin:TrainColMax, :);
        Rs_pre = permute(Rs_pre, [3, 1, 2]);
        Rs_pre = reshape(Rs_pre, [Ncoil, (TrainRowMax-TrainRowMin+1)*(TrainColMax-TrainColMin+1)]);
        Rs = Rs_pre*Rs_pre';
        
        % eigenvalue decomposition
        P = Rn_inv*Rs;
%         P = Rn\Rs;
        [V, D] = eig(P);
        % vector of weights is the eigen-vector of P that corresponds to
        % the largest eigen value
        [~, eigInd] = max(diag(D));
        wt = V(:,eigInd);
        % scale wt to have approximately uniform noise variance
        wt = wt / (wt'*Rn_inv*wt);
        % wt = wt / sqrt(wt'*Rn_inv*wt);
%         wt = wt / sqrt(wt'*Rn\wt);

        % phase correction based on coil with max intensity
        wt = wt.*exp(-sqrt(-1)*argument(wt(maxCoil)));

        wtCoarse(floor((iRow-1)/blocksize)+1, floor((iCol-1)/blocksize)+1,:) = wt;

        if ~bInterp % apply the weights to the whole block
            % need to use the conjugate of the weights, based on Walsh's
            % paper
            BlockRowMin = max(iRow-floor(blocksize/2), 1);
            BlockRowMax = min(iRow-floor(blocksize/2)+blocksize-1, Nrow);
            BlockColMin = max(iCol-floor(blocksize/2), 1);
            BlockColMax = min(iCol-floor(blocksize/2)+blocksize-1, Ncol);
            ImBlock = Im_coil(BlockRowMin:BlockRowMax, BlockColMin:BlockColMax, :);
            ImBlock = permute(ImBlock, [3, 1, 2]);
            ImBlock = reshape(ImBlock, [Ncoil, (BlockRowMax-BlockRowMin+1)*(BlockColMax-BlockColMin+1)]);
            Im(BlockRowMin:BlockRowMax, BlockColMin:BlockColMax) = reshape(wt'*ImBlock, (BlockRowMax-BlockRowMin+1), (BlockColMax-BlockColMin+1));
        end
    end
end

wtFull = zeros(Nrow, Ncol, Ncoil);
if bInterp % interpolate weights to full image matrix, separately for magnitude and phase to avoid phase discontinuity
    for iCoil = 1:Ncoil
        wtFull(:,:,iCoil) = imresize(abs(wtCoarse(:,:,iCoil)),[Nrow, Ncol], 'bilinear') .* exp(sqrt(-1)*imresize(argument(wtCoarse(:,:,iCoil)),[Nrow, Ncol], 'nearest'));
    end
    % need to use the conjugate of the weights, based on Walsh's paper
    Im = sum(Im_coil .* conj(wtFull), 3);
end

