function [F, HalfMaximum] = gaussian_fwhm2d(SIGMA)

% [F] = gaussian_fwhm2d(SIGMA)
%
%	DESCRIPTION
%	Returns a convolution kernel for the 2D covariance matrix SIGMA. SIGMA
%	must be a 2x2 positive-definite covariance matrix. The values of the
%	covariance matrix are cut off at around the half-maximum.

if(~isequal(size(SIGMA), [2, 2]) || ~isreal(SIGMA))
	error('SIGMA must be a 2x2 real matrix');
end

if(SIGMA(1, 2) ~= SIGMA(2, 1))
	error('SIGMA must be symmetric');
end

DetSIGMA = SIGMA(1, 1) * SIGMA(2, 2) - SIGMA(1, 2) * SIGMA(2, 1);

if(DetSIGMA <= 0)
	error('SIGMA must be positive-definite');
end

% PrecisionMatrix = inv(SIGMA)
PrecisionMatrix = [SIGMA(2, 2), -SIGMA(1, 2); -SIGMA(2, 1), SIGMA(1, 1)] ./ DetSIGMA;

XWidth = ceil(abs(SIGMA(1, 1)) / 3);
YWidth = ceil(abs(SIGMA(2, 2)) / 3);
xx = -XWidth:XWidth;
yy = -YWidth:YWidth;

[X, Y] = meshgrid(xx, yy);

XY = [X(:)'; Y(:)'];
clear X Y;

Maximum = 1 ./ ((2 * pi) .* sqrt(DetSIGMA));
HalfMaximum = Maximum / 2;
QuadForm = -0.5 * (XY' * PrecisionMatrix);

F = sum(QuadForm .* XY', 2);
F = reshape(F, length(yy), length(xx));
F = exp(F) .* Maximum;