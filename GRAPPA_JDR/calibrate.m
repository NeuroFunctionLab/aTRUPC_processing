function [kernel,rawkernel] = calibrate(AtA, kSize, nCoil, coil, lambda, sampling) %#codegen
% Dengrong Jiang, July 2018, Johns Hopkins BME
% Modified from calibrate.m in Michael Lustig's SPIRiT_v0.3 package
% Enable 3D GRAPPA kernel calibration

if numel(sampling) == prod(kSize) % if the sampling pattern does not include all coil channels
	sampling = repmat(sampling(:), [nCoil, 1]);
end

% index of the center voxel inside a kernel
if length(kSize) == 2 % 2D GRAPPA
    idxY = (coil-1)*prod(kSize) + sub2ind(kSize, (kSize(1)+1)/2, (kSize(2)+1)/2);
elseif length(kSize) == 3 % 3D GRAPPA
    idxY = (coil-1)*prod(kSize) + sub2ind(kSize, (kSize(1)+1)/2, (kSize(2)+1)/2, (kSize(3)+1)/2);
end
sampling(idxY) = 0;
idxA = find(sampling);

Aty = AtA(:,idxY); Aty = Aty(idxA);
AtA = AtA(idxA,:); AtA =  AtA(:,idxA);

kernel = sampling*0;

lambda = norm(AtA,'fro')/size(AtA,1)*lambda;

rawkernel = (AtA + eye(size(AtA))*lambda)\Aty;
kernel(idxA) = rawkernel; 














