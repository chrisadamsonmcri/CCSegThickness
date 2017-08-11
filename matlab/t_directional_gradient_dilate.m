clear;

[ProbNII, ProbAVW] = load_nii('all_cc_prob.nii.gz');

%H = fspecial('prewitt');

% H = [...
% 	repmat(2, 1, 5); ...
% 	repmat(1, 1, 5); ...
% 	repmat(0, 1, 5); ...
% 	repmat(-1, 1, 5); ...
% 	repmat(-2, 1, 5); ...
% 	];
% H = fspecial('sobel');

F = gaussian_filter_max_1d(3);
ProbAVW = imfilter(ProbAVW, F(:), 'same', 'conv', 'replicate');
ProbAVW = imfilter(ProbAVW, F(:)', 'same', 'conv', 'replicate');

%GX = imfilter(ProbAVW, H', 'same', 'conv', 'replicate');
%GY = imfilter(ProbAVW, H, 'same', 'conv', 'replicate');

BW = ProbAVW < 0.01;

InnerB = bwboundaries(~BW, 4, 'noholes');
OuterB = bwboundaries(imdilate(~BW, ones(3, 3)), 4, 'noholes');

InnerB = InnerB{1};
OuterB = OuterB{1};

InnerBIDX = sub2ind(size(BW), InnerB(:, 1), InnerB(:, 2));
OuterBIDX = sub2ind(size(BW), OuterB(:, 1), OuterB(:, 2));

InnerBSmoothed = imfilter(InnerB, [1, 1, 1, 1, 1, 1 ,1]' / 7, 'same', 'conv', 'circular');
OuterBSmoothed = imfilter(OuterB, [1, 1, 1, 1, 1, 1, 1]' / 7, 'same', 'conv', 'circular');

InnerBGradient = imfilter(InnerBSmoothed, [1, 0, -1]', 'same', 'conv', 'circular');
OuterBGradient = imfilter(OuterBSmoothed, [1, 0, -1]', 'same', 'conv', 'circular');

% clf;
% hold off;
% imshow(~BW, []);
% hold on;
% plot(InnerB(:, 2), InnerB(:, 1), 'r');
% plot(OuterB(:, 2), OuterB(:, 1), 'b');

%hold on;
%plot(BSmoothed(:, 2), BSmoothed(:, 1), 'b')
%BD = conv2([B(end, :); B; B(1, :)], [1 0 -1]', 'valid');

[InnerD, InnerL] = bwdist(BW, 'euclidean');
[OuterD, OuterL] = bwdist(~BW, 'euclidean');

GX = zeros(size(BW));
GY = zeros(size(BW));

[TF, LOC] = ismember(OuterL, InnerBIDX);
GX(TF) = InnerBGradient(LOC(TF), 2);
GY(TF) = InnerBGradient(LOC(TF), 1);
[TF, LOC] = ismember(InnerL, OuterBIDX);
GX(TF) = OuterBGradient(LOC(TF), 2);
GY(TF) = OuterBGradient(LOC(TF), 1);
clf;
% subplot 121;
% imshow(GX, []);
% subplot 122;
% imshow(GY, []);

%ANGLE = atan2(GY, GX);

% imshow(~BW, []);
% hold on;
% plot(InnerB(:, 2), InnerB(:, 1), 'r');
% plot(OuterB(:, 2), OuterB(:, 1), 'b');

% GXUP = imfilter(ProbAVW, [1 -1 0], 'same', 'conv', 'replicate');
% GXDOWN = imfilter(ProbAVW, [0 1 -1], 'same', 'conv', 'replicate');
% GYUP = imfilter(ProbAVW, [1 -1 0]', 'same', 'conv', 'replicate');
% GYDOWN = imfilter(ProbAVW, [0 1 -1]', 'same', 'conv', 'replicate');
% 
% D = abs(GXUP) < abs(GXDOWN);
% GX = zeros(size(ProbAVW));
% GX(D) = GXUP(D);
% GX(~D) = GXDOWN(~D);
% D = abs(GYUP) < abs(GYDOWN);
% GY = zeros(size(ProbAVW));
% GY(D) = GYUP(D);
% GY(~D) = GYDOWN(~D);

%GX = min(abs(GXUP), abs(GXDOWN));
%GY = imfilter(ProbAVW, H, 'same', 'conv', 'replicate');

%[GGX, GGY] = radial_gradient_dilate_fwhm(ProbAVW > 0.5, double(ProbAVW));
% subplot 221;
% imshow(GX, []);
% subplot 222;
% imshow(GY, []);
% 
% subplot 223;
% imshow(GGX, []);
% subplot 224;

% subplot 121;
% imshow(GGY, []); colorbar;
% subplot 122;
% imshow(atan2(GY, GX), []);

% R = regionprops(ProbAVW > 0, 'BoundingBox');
% GY = imcrop(GY, R.BoundingBox);
% GX = imcrop(GX, R.BoundingBox);
% P = imcrop(ProbAVW, R.BoundingBox);
% 
ANGLE = atan2(GY, GX);
ANGLE = mod(ANGLE, 2 * pi);
ANGLE = mod(ANGLE + pi, pi);
%D = ANGLE > pi
%ANGLE = medfilt2(ANGLE, [7, 7]);
%ANGLE = imdilate(ANGLE, ones([3, 3]));
%ANGLE = mod(ANGLE, pi);
%ANGLE(ANGLE < 0) = pi - abs(ANGLE(ANGLE < 0));
V = exp(1i * ANGLE);% .* sqrt(GX .* GX + GY .* GY);
[X, Y] = meshgrid(1:size(V, 2), 1:size(V, 1));
% 
VR = real(V);
VI = imag(V);
% %
% %VR = medfilt2(VR, [5, 5]);
% %VI = medfilt2(VI, [5, 5]);
% clf;
% 
% quiver(X, Y, VR .* (ProbAVW > 0), VI .* (ProbAVW > 0), 0);
% axis ij;
% axis equal ij;

save all_cc_dilation_angle ANGLE;

T = radial_gradient_dilate_fwhm(~BW, ANGLE);

clf;
imshow(T, []);
return;
% ANGLE = atan2(GY, GX) + pi / 2;
% ANGLE = mod(ANGLE, 2 * pi);
% 
% %ANGLE(ANGLE < 0) = pi - abs(ANGLE(ANGLE < 0));
% V = exp(1i * ANGLE);% .* sqrt(GX .* GX + GY .* GY);
% [X, Y] = meshgrid(1:size(V, 2), 1:size(V, 1));
% 
% VR = real(V);
% VI = imag(V);
% 
% subplot 212;
% quiver(X, Y, VR, VI, 0);
% axis equal;
% return;
% subplot 221;
% imshow(GX, []);
% subplot 222;
% imshow(GY, []);
%%

% [X, Y] = getpts;
% 
% I = round(Y);
% J = round(X);

%%

Angle = 45 * pi / 180;

Angle = mod(Angle, 2 * pi);

R = [25, 25] + 25 * abs([cos(Angle), sin(Angle)]);

%[GX(I, J), GY(I, J)]

AngleWeighting = 0.9 * cos(2 * (Angle + 45 * pi / 180));

SQ = sqrt(R(1) * R(2)) * AngleWeighting;

SIGMA = [R(1), SQ; SQ, R(2)];

[F, HalfMaximum] = gaussian_fwhm2d(SIGMA);
%F2 = radial_gradient_dilate_fwhm(SIGMA);

T = padarray(F, floor(50 - size(F) / 2), 0, 'both');
subplot 221;
imshow(T, []);
T = padarray(F > HalfMaximum, floor(50 - size(F) / 2), 0, 'both');
subplot 222;
imshow(T, []);
% T = padarray(F2, floor(50 - size(F2) / 2), 0, 'both');
% subplot 223;
% imshow(T, []);
% T = padarray(F2 > HalfMaximum, floor(50 - size(F2) / 2), 0, 'both');
% subplot 224;
% imshow(T, []);
