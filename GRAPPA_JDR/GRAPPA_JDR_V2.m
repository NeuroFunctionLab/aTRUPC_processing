function filledkData = GRAPPA_JDR_V2(kData,kCalib,kernelSize,lambda, Rfactor)
% Dengrong Jiang, 09-22-2015, Johns Hopkins BME
% Modified from GRAPPA.m written by Michael Lustig 2008
[Nky,Nkx,Ncoil] = size(kData);

% decide which direction is the foldover direction, transpose kData if it's
% not along Ky
isKx = sum(abs(kData(round(Nky/2),:,1))>0) < Nkx; % is foldover along kx?
if isKx
    kData = permute(kData,[2,1,3]); % kData must be transposed back before return
    kernelSize = fliplr(kernelSize); % also need to switch the dimensions of kernelSize
    kCalib = permute(kCalib,[2,1,3]);
end

[AtA] = corrMatrix(kCalib,kernelSize); % build coil correlation matrix
	
filledkData = PAT(kData, AtA,kernelSize, lambda, Rfactor); % reconstruct single coil image

if isKx
    filledkData = permute(filledkData,[2,1,3]); % kData must be transposed back before return
end

function [filled_kData] = PAT(kData, AtA, kernelSize, lambda, Rfactor)
% Customized for my own implementation of GRAPPA subsampling
% Dengrong Jiang, JHU BME, 2015-09-15
[Nky,Nkx,Ncoil] = size(kData);

% % decide which direction is the foldover direction, transpose kData if it's
% % not along Ky
% isKx = sum(abs(kData(round(Nky/2),:,iCoil))>0) < Nkx; % is foldover along kx?
% if isKx
%     kData = permute(kData,[2,1,3]); % kData must be transposed back before return
%     [Nky,Nkx,Ncoil] = size(kData); % get the size of current kData
%     kernelSize = fliplr(kernelSize); % also need to switch the dimensions of kernelSize
% end
iCoil = 1;
LineIndicator = abs(kData(:,round(Nkx/2),iCoil))>0;
filled_kData = kData; % copy the sampled points into the k-space to fill
RefUpperEdge = find((LineIndicator(2:end)+LineIndicator(1:end-1))==2,1,'first');
RefLowerEdge = find((LineIndicator(2:end)+LineIndicator(1:end-1))==2,1,'last')+1;
kernels = zeros(Ncoil,kernelSize(1)*kernelSize(2)*Ncoil);
%% Major part (not including the margins and the lines very close to ref region)    
for iRfactor = 1:Rfactor-1
    % calculate kernel
    centerky = RefLowerEdge+Rfactor+iRfactor; % the center of the kernel
    centerkx = round(Nkx/2);
    tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
    pattern = abs(tmp)>0;
    for iCoil = 1:Ncoil
        kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
        kernels(iCoil,:) = transpose(kernel_iCoil(:));
    end
    % extract ky and kx involved
    kyLinToFill = [fliplr((RefUpperEdge-2*Rfactor+iRfactor):(-Rfactor):1+(kernelSize(1)-1)/2),...
        (RefLowerEdge+Rfactor+iRfactor):Rfactor:Nky-(kernelSize(1)-1)/2];
    kxToFill = [((kernelSize(2)-1)/2+1) : Nkx-(kernelSize(2)-1)/2];
    %% construct index matrix for the points used to interpolate missing
    % points
    % part1: represents the index increment along a ky line
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    part1 = repmat(0:Nky:Nky*(length(kxToFill)-1),[kernelSize(1)*kernelSize(2)*Ncoil, length(kyLinToFill)]);
    % part2: represents the index increment through ky lines
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    part2 = repmat(kyLinToFill - kyLinToFill(1),[length(kxToFill),1]);
    part2 = repmat(transpose(part2(:)),[kernelSize(1)*kernelSize(2)*Ncoil, 1]);
    % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
    part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
    part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
    part3_2 = part3_2(:);
    part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
    part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
    part3_4 = part3_4(:);
    part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
    % Final index matrix for the points used to interpolate missing
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part2 + part1;
    %% Fill missing points and assign them to proper position in filled_kData
    filled_points = kernels * kData(IndMatrix);
    filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
    filled_kData(kyLinToFill,kxToFill,:) = permute(filled_points,[2,1,3]);
end

%% Right above RefUpperEdge
for iRfactor = 1:Rfactor-1
    % calculate kernel
    centerky = RefUpperEdge-iRfactor; % the center of the kernel
    centerkx = round(Nkx/2);
    tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
    pattern = abs(tmp)>0;
    for iCoil = 1:Ncoil
        kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
        kernels(iCoil,:) = transpose(kernel_iCoil(:));
    end
    % extract ky and kx involved
    kyLinToFill = centerky;
    kxToFill = [((kernelSize(2)-1)/2+1) : Nkx-(kernelSize(2)-1)/2];
    %% construct index matrix for the points used to interpolate missing
    % points
    % part1: represents the index increment along a ky line
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    part1 = repmat(0:Nky:Nky*(length(kxToFill)-1),[kernelSize(1)*kernelSize(2)*Ncoil, length(kyLinToFill)]);
    % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
    part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
    part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
    part3_2 = part3_2(:);
    part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
    part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
    part3_4 = part3_4(:);
    part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
    % Final index matrix for the points used to interpolate missing
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part1;
    %% Fill missing points and assign them to proper position in filled_kData
    filled_points = kernels * kData(IndMatrix);
    filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
    filled_kData(kyLinToFill,kxToFill,:) = permute(filled_points,[2,1,3]);
end

%% Right below RefLowerEdge
for iRfactor = 1:Rfactor-1
    % calculate kernel
    centerky = RefLowerEdge+iRfactor; % the center of the kernel
    centerkx = round(Nkx/2);
    tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
    pattern = abs(tmp)>0;
    for iCoil = 1:Ncoil
        kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
        kernels(iCoil,:) = transpose(kernel_iCoil(:));
    end
    % extract ky and kx involved
    kyLinToFill = centerky;
    kxToFill = [((kernelSize(2)-1)/2+1) : Nkx-(kernelSize(2)-1)/2];
    %% construct index matrix for the points used to interpolate missing
    % points
    % part1: represents the index increment along a ky line
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    part1 = repmat(0:Nky:Nky*(length(kxToFill)-1),[kernelSize(1)*kernelSize(2)*Ncoil, length(kyLinToFill)]);
    % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
    part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
    part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
    part3_2 = part3_2(:);
    part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
    part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
    part3_4 = part3_4(:);
    part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
    % Final index matrix for the points used to interpolate missing
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part1;
    %% Fill missing points and assign them to proper position in filled_kData
    filled_points = kernels * kData(IndMatrix);
    filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
    filled_kData(kyLinToFill,kxToFill,:) = permute(filled_points,[2,1,3]);
end

%% Zero padding the surrounding of k-space
kData = zpad(kData,[Nky+kernelSize(1)-1, Nkx+kernelSize(2)-1,Ncoil]);
[Nky, Nkx, Ncoil] = size(kData);

%% Upper margin
for ikernelSize = 1:(kernelSize(1)-1)/2
    % calculate kernel
    centerky = kernelSize(1) - ikernelSize; % the center of the kernel
    centerkx = round(Nkx/2);
    if sum(abs(kData(centerky,centerkx,1))) > 0 % if the line is sampled, no GRAPPA is needed
        continue;
    end
    tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
    pattern = abs(tmp)>0;
    for iCoil = 1:Ncoil
        kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
        kernels(iCoil,:) = transpose(kernel_iCoil(:));
    end
    % extract ky and kx involved
    kyLinToFill = centerky;
    kxToFill = [kernelSize(2):(Nkx-kernelSize(2)+1)]; % Note the size of kData has changed
    %% construct index matrix for the points used to interpolate missing
    % points
    % part1: represents the index increment along a ky line
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    part1 = repmat(0:Nky:Nky*(length(kxToFill)-1),[kernelSize(1)*kernelSize(2)*Ncoil, length(kyLinToFill)]);
    % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
    part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
    part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
    part3_2 = part3_2(:);
    part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
    part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
    part3_4 = part3_4(:);
    part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
    % Final index matrix for the points used to interpolate missing
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part1;
    %% Fill missing points and assign them to proper position in filled_kData
    filled_points = kernels * kData(IndMatrix);
    filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
    orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
    orikxToFill = kxToFill - (kernelSize(2)-1)/2;
    filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
end

%% Lower margin
for ikernelSize = 1:(kernelSize(1)-1)/2
    % calculate kernel
    centerky = size(kData,1)-kernelSize(1)+1 + ikernelSize; % the center of the kernel
    centerkx = round(Nkx/2);
    if sum(abs(kData(centerky,centerkx,1))) > 0 % if the line is sampled, no GRAPPA is needed
        continue;
    end
    tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
    pattern = abs(tmp)>0;
    for iCoil = 1:Ncoil
        kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
        kernels(iCoil,:) = transpose(kernel_iCoil(:));
    end
    % extract ky and kx involved
    kyLinToFill = centerky;
    kxToFill = [kernelSize(2):(Nkx-kernelSize(2)+1)]; % Note the size of kData has changed
    %% construct index matrix for the points used to interpolate missing
    % points
    % part1: represents the index increment along a ky line
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    part1 = repmat(0:Nky:Nky*(length(kxToFill)-1),[kernelSize(1)*kernelSize(2)*Ncoil, length(kyLinToFill)]);
    % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
    part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
    part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
    part3_2 = part3_2(:);
    part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
    part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
    part3_4 = part3_4(:);
    part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
    % Final index matrix for the points used to interpolate missing
    % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
    IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part1;
    %% Fill missing points and assign them to proper position in filled_kData
    filled_points = kernels * kData(IndMatrix);
    filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
    orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
    orikxToFill = kxToFill - (kernelSize(2)-1)/2;
    filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
end

%% Left margin of major part
for ikernelSize = 1:(kernelSize(2)-1)/2
    for iRfactor = 1:Rfactor-1    
        % calculate kernel
        centerky = (kernelSize(1)-1)/2+RefLowerEdge+Rfactor+iRfactor; % Note the size change of kData
        centerkx = (kernelSize(1)-1)/2 + ikernelSize;
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        % extract ky and kx involved
        kyLinToFill = [fliplr(((kernelSize(1)-1)/2+RefUpperEdge-2*Rfactor+iRfactor):(-Rfactor):kernelSize(1)),...
            ((kernelSize(1)-1)/2+RefLowerEdge+Rfactor+iRfactor):Rfactor:(Nky-kernelSize(1)+1)]; % Note the size change of kData
        kxToFill = centerkx;
        %% construct index matrix for the points used to interpolate missing
        % points
        % part2: represents the index increment through ky lines
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        part2 = repmat(kyLinToFill - kyLinToFill(1),[length(kxToFill),1]);
        part2 = repmat(transpose(part2(:)),[kernelSize(1)*kernelSize(2)*Ncoil, 1]);
        % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
        part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
        part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
        part3_2 = part3_2(:);
        part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
        part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
        part3_4 = part3_4(:);
        part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
        % Final index matrix for the points used to interpolate missing
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part2;
        %% Fill missing points and assign them to proper position in filled_kData
        filled_points = kernels * kData(IndMatrix);
        filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
        orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = kxToFill - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
    end
end

%% Right margin of major part
for iRfactor = 1:Rfactor-1
    for ikernelSize = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = (kernelSize(1)-1)/2+RefLowerEdge+Rfactor+iRfactor; % the center of the kernel
        centerkx = size(kData,2)-kernelSize(1)+1 + ikernelSize;
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        % extract ky and kx involved
        kyLinToFill = [fliplr(((kernelSize(1)-1)/2+RefUpperEdge-2*Rfactor+iRfactor):(-Rfactor):kernelSize(1)),...
            ((kernelSize(1)-1)/2+RefLowerEdge+Rfactor+iRfactor):Rfactor:(Nky-kernelSize(1)+1)]; % Note the size change of kData
        kxToFill = centerkx;
        %% construct index matrix for the points used to interpolate missing
        % points
        % part2: represents the index increment through ky lines
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        part2 = repmat(kyLinToFill - kyLinToFill(1),[length(kxToFill),1]);
        part2 = repmat(transpose(part2(:)),[kernelSize(1)*kernelSize(2)*Ncoil, 1]);
        % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
        part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
        part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
        part3_2 = part3_2(:);
        part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
        part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
        part3_4 = part3_4(:);
        part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
        % Final index matrix for the points used to interpolate missing
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part2;
        %% Fill missing points and assign them to proper position in filled_kData
        filled_points = kernels * kData(IndMatrix);
        filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
        orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = kxToFill - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
    end
end

%% Left margin right above RefUpperEdge
for iRfactor = 1:Rfactor-1
    for ikernelSize = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = (kernelSize(1)-1)/2+RefUpperEdge-iRfactor; % Note the size change of kData
        centerkx = kernelSize(1) - ikernelSize;
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        % extract ky and kx involved
        kyLinToFill = centerky; 
        kxToFill = centerkx;
        %% construct index matrix for the points used to interpolate missing
        % points
        % part2: represents the index increment through ky lines
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        part2 = repmat(kyLinToFill - kyLinToFill(1),[length(kxToFill),1]);
        part2 = repmat(transpose(part2(:)),[kernelSize(1)*kernelSize(2)*Ncoil, 1]);
        % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
        part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
        part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
        part3_2 = part3_2(:);
        part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
        part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
        part3_4 = part3_4(:);
        part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
        % Final index matrix for the points used to interpolate missing
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part2;
        %% Fill missing points and assign them to proper position in filled_kData
        filled_points = kernels * kData(IndMatrix);
        filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
        orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = kxToFill - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
    end
end

%% Right margin right above RefUpperEdge
for iRfactor = 1:Rfactor-1
    for ikernelSize = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = (kernelSize(1)-1)/2+RefUpperEdge-iRfactor; % Note the size change of kData
        centerkx = size(kData,2)-kernelSize(1)+1 + ikernelSize;
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        % extract ky and kx involved
        kyLinToFill = centerky; 
        kxToFill = centerkx;
        %% construct index matrix for the points used to interpolate missing
        % points
        % part2: represents the index increment through ky lines
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        part2 = repmat(kyLinToFill - kyLinToFill(1),[length(kxToFill),1]);
        part2 = repmat(transpose(part2(:)),[kernelSize(1)*kernelSize(2)*Ncoil, 1]);
        % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
        part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
        part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
        part3_2 = part3_2(:);
        part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
        part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
        part3_4 = part3_4(:);
        part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
        % Final index matrix for the points used to interpolate missing
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part2;
        %% Fill missing points and assign them to proper position in filled_kData
        filled_points = kernels * kData(IndMatrix);
        filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
        orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = kxToFill - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
    end
end

%% Left margin right below RefLowerEdge
for iRfactor = 1:Rfactor-1
    for ikernelSize = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = (kernelSize(1)-1)/2+RefLowerEdge+iRfactor; % Note the size change of kData
        centerkx = kernelSize(1) - ikernelSize;
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        % extract ky and kx involved
        kyLinToFill = centerky; 
        kxToFill = centerkx;
        %% construct index matrix for the points used to interpolate missing
        % points
        % part2: represents the index increment through ky lines
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        part2 = repmat(kyLinToFill - kyLinToFill(1),[length(kxToFill),1]);
        part2 = repmat(transpose(part2(:)),[kernelSize(1)*kernelSize(2)*Ncoil, 1]);
        % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
        part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
        part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
        part3_2 = part3_2(:);
        part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
        part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
        part3_4 = part3_4(:);
        part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
        % Final index matrix for the points used to interpolate missing
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part2;
        %% Fill missing points and assign them to proper position in filled_kData
        filled_points = kernels * kData(IndMatrix);
        filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
        orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = kxToFill - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
    end
end

%% Right margin right below RefLowerEdge
for iRfactor = 1:Rfactor-1
    for ikernelSize = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = (kernelSize(1)-1)/2+RefLowerEdge+iRfactor; % the center of the kernel
        centerkx = size(kData,2)-kernelSize(1)+1 + ikernelSize;
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        % extract ky and kx involved
        kyLinToFill = centerky; 
        kxToFill = centerkx;
        %% construct index matrix for the points used to interpolate missing
        % points
        % part2: represents the index increment through ky lines
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        part2 = repmat(kyLinToFill - kyLinToFill(1),[length(kxToFill),1]);
        part2 = repmat(transpose(part2(:)),[kernelSize(1)*kernelSize(2)*Ncoil, 1]);
        % part3: indices of the kernel centered at [kyLinToFill(1),kxToFill(1)]
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X 1
        part3_1 = [kyLinToFill(1)-(kernelSize(1)-1)/2:kyLinToFill(1)+(kernelSize(1)-1)/2].'+(kxToFill(1)-1-(kernelSize(2)-1)/2)*Nky; % first column
        part3_2 = repmat([0:Nky:(Nky*(kernelSize(2)-1))],[kernelSize(1),1]); % index increment all columns
        part3_2 = part3_2(:);
        part3_3 = repmat(part3_1,[kernelSize(2),1]) + part3_2; % index of kernel in Coil 1
        part3_4 = repmat([0:Nky*Nkx:Nky*Nkx*(Ncoil-1)],[kernelSize(1)*kernelSize(2),1]);
        part3_4 = part3_4(:);
        part3 = repmat(part3_3,[Ncoil,1]) + part3_4;
        % Final index matrix for the points used to interpolate missing
        % Size: (kernelSize(1)*kernelSize(2)*Ncoil) X (length(kxToFill)*length(kyLinToFill))
        IndMatrix = repmat(part3,[1,length(kxToFill)*length(kyLinToFill)]) + part2;
        %% Fill missing points and assign them to proper position in filled_kData
        filled_points = kernels * kData(IndMatrix);
        filled_points = reshape(transpose(filled_points), [length(kxToFill),length(kyLinToFill),Ncoil]);
        orikyLinToFill = kyLinToFill - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = kxToFill - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = permute(filled_points,[2,1,3]);
    end
end

%% Left upper corner
for ikernelSizeKy = 1:(kernelSize(1)-1)/2
    for ikernelSizeKx = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = kernelSize(1) - ikernelSizeKy; % the center of the kernel
        centerkx = kernelSize(2) - ikernelSizeKx;
        if sum(abs(kData(centerky,centerkx,1))) > 0 % if the line is sampled, no GRAPPA is needed
            continue;
        end        
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        filled_points = kernels * tmp(:);
        filled_points = reshape(filled_points, [length(kyLinToFill),length(kxToFill),Ncoil]);
        orikyLinToFill = centerky - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = centerkx - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = filled_points;
    end
end

%% Right upper corner
for ikernelSizeKy = 1:(kernelSize(1)-1)/2
    for ikernelSizeKx = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = kernelSize(1) - ikernelSizeKy; % the center of the kernel
        centerkx = Nkx - kernelSize(2) + 1 + ikernelSizeKx;
        if sum(abs(kData(centerky,centerkx,1))) > 0 % if the line is sampled, no GRAPPA is needed
            continue;
        end        
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        filled_points = kernels * tmp(:);
        filled_points = reshape(filled_points, [length(kyLinToFill),length(kxToFill),Ncoil]);
        orikyLinToFill = centerky - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = centerkx - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = filled_points;
    end
end

%% Left lower corner
for ikernelSizeKy = 1:(kernelSize(1)-1)/2
    for ikernelSizeKx = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = Nky - kernelSize(1) + 1 + ikernelSizeKy; % the center of the kernel
        centerkx = kernelSize(2) - ikernelSizeKx;
        if sum(abs(kData(centerky,centerkx,1))) > 0 % if the line is sampled, no GRAPPA is needed
            continue;
        end          
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        filled_points = kernels * tmp(:);
        filled_points = reshape(filled_points, [length(kyLinToFill),length(kxToFill),Ncoil]);
        orikyLinToFill = centerky - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = centerkx - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = filled_points;
    end
end

%% Right lower corner
for ikernelSizeKy = 1:(kernelSize(1)-1)/2
    for ikernelSizeKx = 1:(kernelSize(2)-1)/2
        % calculate kernel
        centerky = Nky - kernelSize(1) + 1 + ikernelSizeKy; % the center of the kernel
        centerkx = Nkx - kernelSize(2) + 1 + ikernelSizeKx;
        if sum(abs(kData(centerky,centerkx,1))) > 0 % if the line is sampled, no GRAPPA is needed
            continue;
        end          
        tmp = kData(centerky-(kernelSize(1)-1)/2:centerky+(kernelSize(1)-1)/2,centerkx-(kernelSize(2)-1)/2:centerkx+(kernelSize(2)-1)/2,:);
        pattern = abs(tmp)>0;
        for iCoil = 1:Ncoil
            kernel_iCoil = calibrate(AtA,kernelSize,Ncoil,iCoil,lambda,pattern);
            kernels(iCoil,:) = transpose(kernel_iCoil(:));
        end
        filled_points = kernels * tmp(:);
        filled_points = reshape(filled_points, [length(kyLinToFill),length(kxToFill),Ncoil]);
        orikyLinToFill = centerky - (kernelSize(1)-1)/2; % Note the size change of kData
        orikxToFill = centerkx - (kernelSize(2)-1)/2;
        filled_kData(orikyLinToFill,orikxToFill,:) = filled_points;
    end
end

% if isKx
%     filled_kData = permute(filled_kData,[2,1,3]); % kData must be transposed back before return
% end
