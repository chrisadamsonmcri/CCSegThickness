function [BWI] = diagonal_dilate(BW, R, AngleWeighting)

% [F] = diagonal_dilate(BW, R)
%
% DESCRIPTION
%	
% PARAMETERS
%
% RETURNS

SQ = sqrt(R(1) * R(2)) * AngleWeighting;

SIGMA = [R(1), SQ; SQ, R(2)];

[F] = gaussian_fwhm2d(SIGMA);

[ma,na] = size(BW);
[mb,nb] = size(F);

szm = max([ma+mb-1,ma,mb]);
szn = max([na+nb-1,na,nb]);

FF = fft2(F, szm, szn);
BWF = fft2(double(BW), szm, szn);

BWI = real(ifft2(FF .* BWF));

%BWI = BWI(size(F, 1):size(F, 1) + size(BW, 1) - 1, size(F, 2):size(F, 2) + size(BW, 2) - 1);

px = ((size(F, 2) - 1) + mod((size(F, 2) - 1), 2)) / 2;
py = ((size(F, 1) - 1) + mod((size(F, 1) - 1), 2)) / 2;

BWI = BWI(py + 1:py + size(BW, 1), px + 1:px + size(BW, 2));