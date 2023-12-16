function a = argument(z)
   re = real(z);
   im = imag(z);
   
   a=atan2(im,re);
% Dengrong June 2020: the codes below will give incorrect outputs if 
% the real part is exactly 0! If the imaginary part is not exactly 0, it 
% will output pi; if the imaginary part is exactly 0, it will output NaN! 
% This would be a problem for masked images. 
% For complex images with no voxel whose real part is exactly zero, the 
% difference between the phase computed by the codes below and that by 
% atan2 is less than 2*eps (the Matlab machine precision).
% Therefore, it is better to just use atan2

    %a=atan(im/re)+(pi/2)*sign(im)*(1-sign(re));
%   a=atan(im./re)+(pi/2)*sign(im).*(ones(size(re))-sign(re));
  %a=atan(im/re)+(pi/2)*sign(im)*(1-sign(re));