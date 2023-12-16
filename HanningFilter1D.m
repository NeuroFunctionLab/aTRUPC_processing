function [flt, FWHM, PSF] = HanningFilter1D(Ntotal, Nflat, Nslope, display_expand)
% This function computes the 1D Hanning Filter
% Input:
% - Ntotal: total number of points
% - Nflat: length of plateau
% - Nslope: length of slope
% - display_expand: zero-padding to a factor of display_expand to compute PSF
% Output:
% - flt: Hanning filter
% - FWHM: FWHM of the PSF
% - PSF: point spread function
% Dengrong Jiang, JHU BME, Feb 2018

if nargin < 4
    display_expand = 500;
end
if Nflat+2*Nslope > Ntotal
    fprintf('Error: Nflat+2*Nslope > Ntotal\n');
end
flt = zeros(1, Ntotal);
midpoint = floor((Ntotal+1)/2);
leftend = midpoint-floor(Nflat/2)+1;
rightend = midpoint-floor(Nflat/2)+Nflat;
flt(leftend:rightend) = 1;
% flt(leftend-Nslope:leftend-1) = 1/2*(1+cos(pi*[Nslope:-1:1]/Nslope));
% flt(rightend+1:rightend+Nslope) = 1/2*(1+cos(pi*[1:Nslope]/Nslope));
flt(leftend-Nslope:leftend-1) = 1/2*(1+cos(pi*[Nslope:-1:1]/(Nslope+1)));
flt(rightend+1:rightend+Nslope) = 1/2*(1+cos(pi*[1:Nslope]/(Nslope+1)));
% compute PSF and FWHM
PSF = fftshift(ifft(ifftshift(zpad(flt, 1, display_expand*length(flt)))));
FWHM = length(find(abs(PSF)>=1/2*max(abs(PSF)))) / display_expand;
