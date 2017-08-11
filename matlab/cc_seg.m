clear;

[FornixProbNII, FornixProbAVW] = load_nii('all_fornix_prob.nii.gz');
[ProbNII, ProbAVW] = load_nii('all_cc_prob.nii.gz');
[TemplateNII, TemplateAVW] = load_nii('all_msp_mean.nii.gz');
[SmoothNII, SmoothAVW] = load_nii('all_msp_smooth_areas.nii.gz');

SmoothAVW = (SmoothAVW > 0);
SmoothAVW = imerode(SmoothAVW, ones(3, 3));
F = gaussian_filter_max_1d(1);

SmoothAVW = double(SmoothAVW);

SmoothAVW = imfilter(SmoothAVW, F(:), 'same', 'conv', 'replicate');
SmoothAVW = imfilter(SmoothAVW, F(:)', 'same', 'conv', 'replicate');
%SmoothAVW
%keyboard;
SmoothAVW = SmoothAVW ./ max(SmoothAVW(:));
SmoothAVW = SmoothAVW .* 3 + 1;
[NII, AVW] = load_nii('0001_acpc.nii.gz');
%keyboard;
%[RegNII, RegAVW] = load_nii('0001_acpc_reg.nii.gz');
TemplatePixdims = TemplateNII.hdr.dime.pixdim(3:4);
NIIPixdims = NII.hdr.dime.pixdim(3:4);

MidsagSliceIDX = size(AVW, 2) / 2;
XI = floor(MidsagSliceIDX);
XFrac = MidsagSliceIDX - XI;

%[Y, Z] = meshgrid(size(AVW, 1):-1:1, size(AVW, 3):-1:1);

%AVWMidSag = interp3(double(AVW), repmat(MidsagSliceIDX, size(Y)), Y, Z);
%clear Y Z;
AVWMidSag = AVW(size(AVW, 1):-1:1, XI, size(AVW, 3):-1:1) .* (1 - XFrac) + ...
	AVW(size(AVW, 1):-1:1, XI + 1, size(AVW, 3):-1:1) .* XFrac;
AVWMidSag = squeeze(AVWMidSag)';
AVWMidSag = double(AVWMidSag);
clear XI XFrac MidsagSliceIDX;

% resample to the space of the template
AVWMidSagxx = (1:size(AVWMidSag, 2)) * NIIPixdims(1);
AVWMidSagxx = AVWMidSagxx - mean(AVWMidSagxx);

AVWMidSagyy = (1:size(AVWMidSag, 1)) * NIIPixdims(2);
AVWMidSagyy = AVWMidSagyy - mean(AVWMidSagyy);

[AVWMidSagX, AVWMidSagY] = meshgrid(AVWMidSagxx, AVWMidSagyy);

Templatexx = (1:size(TemplateAVW, 2)) * TemplatePixdims(1);
Templatexx = Templatexx - mean(Templatexx);

Templateyy = (1:size(TemplateAVW, 1)) * TemplatePixdims(2);
Templateyy = Templateyy - mean(Templateyy);

[TemplateX, TemplateY] = meshgrid(Templatexx, Templateyy);

%R = [163.7248  168.6860  201.4264  125.0233];
%R = [187.3418  195.7478  160.4776   65.7194];
R = [184.1436  189.6490  165.1613   92.8049];
MatchedFilter = imcrop(TemplateAVW, R);
MatchedFilterProb = imcrop(ProbAVW, R);
MatchedFilterFornixProb = imcrop(FornixProbAVW, R);
MatchedFilterSmoothAVW = imcrop(SmoothAVW, R);

LeftLimit = max(min(Templatexx), min(AVWMidSagxx));
UpperLimit = max(min(Templateyy), min(AVWMidSagyy));
RightLimit = min(max(Templatexx), max(AVWMidSagxx));
LowerLimit = min(max(Templateyy), max(AVWMidSagyy));

[ResampleX, ResampleY] = meshgrid(LeftLimit:TemplatePixdims(1):RightLimit, UpperLimit:TemplatePixdims(2):LowerLimit);

ResampledAVW = interp2(AVWMidSagX, AVWMidSagY, AVWMidSag, ResampleX, ResampleY);

cc = normxcorr2(MatchedFilter, ResampledAVW);
[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));
corr_offset = [ (ypeak-size(MatchedFilter,1)) (xpeak-size(MatchedFilter,2)) ];
%corr_offset = round(corr_offset + R([2 1]));

%rect_offset = [size(ResampledAVW, 1) - size(MatchedFilter, 1), size(ResampledAVW, 2) - size(MatchedFilter, 2)];

total_offset = round(corr_offset);
E = edge(MatchedFilter, 'canny', [0.1 0.4]);
%E(:) = 1;
EIMG = repmat(reshape([0, 0, 1], [1, 1, 3]), [size(ResampledAVW), 1]);
AlphaData = zeros(size(ResampledAVW));
%AlphaData = E;
%F = ResampledAVW;
AlphaData(total_offset(1):total_offset(1) + size(MatchedFilter, 1) - 1, total_offset(2):total_offset(2) + size(MatchedFilter, 2) - 1) = E;

% do histogram matching
InitialResampledAVW = ResampledAVW(total_offset(1):total_offset(1) + size(MatchedFilter, 1) - 1, total_offset(2):total_offset(2) + size(MatchedFilter, 2) - 1);

InitialResampledAVWBins = linspace(min(InitialResampledAVW(:)), max(InitialResampledAVW(:)), 998);
S = InitialResampledAVWBins(2) - InitialResampledAVWBins(1);
InitialResampledAVWBins = [InitialResampledAVWBins(1) - S, InitialResampledAVWBins, InitialResampledAVWBins(end) + S];
clear S;

% get the CDFs of the histogram of the target image
[InitialResampledAVWHist, InitialResampledAVWBins] = hist(InitialResampledAVW(:), InitialResampledAVWBins);
InitialResampledAVWHist = InitialResampledAVWHist ./ sum(InitialResampledAVWHist);
InitialResampledAVWHistCDF = cumsum(InitialResampledAVWHist);
% get the CDFs of the histogram of the matched filter

MatchedFilterBins = linspace(min(MatchedFilter(:)), max(MatchedFilter(:)), 998);
S = MatchedFilterBins(2) - MatchedFilterBins(1);
MatchedFilterBins = [MatchedFilterBins(1) - S, MatchedFilterBins, MatchedFilterBins(end) + S];
clear S;
[MatchedFilterHist, MatchedFilterBins] = hist(MatchedFilter(:), MatchedFilterBins);
MatchedFilterHist = MatchedFilterHist ./ sum(MatchedFilterHist);
MatchedFilterHistCDF = cumsum(MatchedFilterHist);

% use unique to create a CDF with no repeated elements
[UMatchedFilterHistCDF, I, J] = unique(MatchedFilterHistCDF);
UMatchedFilterHistCDFBins = MatchedFilterBins(I);

[UInitialResampledAVWHistCDF, I, J] = unique(InitialResampledAVWHistCDF);
UInitialResampledAVWHistCDFBins = InitialResampledAVWBins(I);
clear I J;

CDFMatchedFilterPixels = interp1(UMatchedFilterHistCDFBins, UMatchedFilterHistCDF, MatchedFilter(:));
MatchedFilterRemapped = interp1(UInitialResampledAVWHistCDF, UInitialResampledAVWHistCDFBins, CDFMatchedFilterPixels);

MatchedFilterRemapped = reshape(MatchedFilterRemapped, size(MatchedFilter));

%plot(InitialResampledAVWBins, InitialResampledAVWHistCDF, MatchedFilterBins, MatchedFilterHistCDF);
% subplot 121;
% imshow(MatchedFilterRemapped, []);
% colorbar;
% subplot 122;
% imshow(InitialResampledAVW, []);
% colorbar;

%return;
disp('Running LK');

%I = (ResampledAVW - min(ResampledAVW(:))) ./ (max(ResampledAVW(:)) - min(ResampledAVW(:)));
%J = (MatchedFilter - min(MatchedFilter(:))) ./ (max(MatchedFilter(:)) - min(MatchedFilter(:)));
Parameters = lk_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, [0 0 0 0 corr_offset(2) corr_offset(1)], 100);

[TemplateXISO, TemplateYISO] = meshgrid(1:size(MatchedFilter, 2), 1:size(MatchedFilter, 1));
%TemplateXISO = TemplateXISO - corr_offset(2);
%TemplateYISO = TemplateYISO - corr_offset(1);
XY = cat(1, TemplateXISO(:)', TemplateYISO(:)', ones(1, numel(MatchedFilter)));
%%
TransformationMatrix = [1 + Parameters(1), Parameters(3), Parameters(5); ...
    Parameters(2), 1 + Parameters(4), Parameters(6); ...
    0, 0, 1];

TXY = TransformationMatrix * XY;

TXY = TXY(1:2, :);

TX = reshape(TXY(1, :), size(MatchedFilter));
TY = reshape(TXY(2, :), size(MatchedFilter));
clear TXY TransformationMatrix;

% find the change in vector along the X direction
XChangeVector = [TX(1, 2) - TX(1, 1), TY(1, 2) - TY(1, 1)];
YChangeVector = [TX(2, 1) - TX(1, 1), TY(2, 1) - TY(1, 1)];

UnitXChangeVector = XChangeVector ./ norm(XChangeVector);
UnitYChangeVector = YChangeVector ./ norm(YChangeVector);	

% find the change in vector along the X direction
% scalar resolute in the direction of XChangeVector

[ResampledAVWX, ResampledAVWY] = meshgrid(1:size(ResampledAVW, 2), 1:size(ResampledAVW, 1));

XDot = (ResampledAVWX - TX(1, 1)) .* UnitXChangeVector(1) + (ResampledAVWY - TY(1, 1)) .* UnitXChangeVector(2);
YDot = (ResampledAVWX - TX(1, 1)) .* UnitYChangeVector(1) + (ResampledAVWY - TY(1, 1)) .* UnitYChangeVector(2);

X = XDot ./ norm(XChangeVector);
Y = YDot ./ norm(YChangeVector);
% interpolate 
%WarpedIMG = interp2(TemplateAVW, TX, TY, 'bilinear');
TIMG = interp2(MatchedFilter, X, Y, 'linear', 0);
PIMG = interp2(MatchedFilterProb, X, Y, 'linear', 0);
FIMG = interp2(MatchedFilterFornixProb, X, Y, 'linear', 0);
SIMG = interp2(MatchedFilterSmoothAVW, X, Y, 'linear', 0);
% E = edge(TIMG, 'canny');
% IMG = (ResampledAVW - min(ResampledAVW(:))) ./ (max(ResampledAVW(:)) - min(ResampledAVW(:)));
% IMG = repmat(IMG, [1, 1, 3]);
% 
% [I, J] = find(E);
% for z = 1:length(I)
% 	IMG(I(z), J(z), :) = [0, 0, 1];
% end
% 
% subplot 121;
% imshow(IMG, []);
% subplot 122;
% imshow(ResampledAVW, []);
% hold on;
% imagesc(EIMG, 'AlphaData', AlphaData);
%WarpedIMG(isnan(WarpedIMG)) = 0;

% find the bounding box of the template warped to the image

R = regionprops(TIMG > 0, 'BoundingBox');

ResampledAVWCropped = imcrop(ResampledAVW, R.BoundingBox);
PIMGCropped = imcrop(PIMG, R.BoundingBox);
FIMGCropped = imcrop(FIMG, R.BoundingBox);
SIMGCropped = imcrop(SIMG, R.BoundingBox);
ResampledAVWCroppedForOtsu = (ResampledAVWCropped - min(ResampledAVWCropped(:))) ./ (max(ResampledAVWCropped(:)) - min(ResampledAVWCropped(:)));
ResampledAVWCroppedForOtsu = round(ResampledAVWCroppedForOtsu * 255);
%%
T = otsu2(ResampledAVWCroppedForOtsu);

OtsuMask = (ResampledAVWCroppedForOtsu >= T(end));

OtsuMask = imopen(OtsuMask, strel('diamond', 2));

L = bwlabel(OtsuMask);
R = regionprops(L, 'Area');

[~, I] = max([R.Area]);

OtsuMask = (L == I);
clear Junk;
%%
% clf;
% subplot 321;
% imshow(ResampledAVWCropped, []);
% subplot 323;
% imshow(ResampledAVWCropped, []);
% hold on;
% [I, J] = find(OtsuMask & (PIMGCropped > 0.5));
% plot(J, I, '*');
% subplot(3, 2, 2);
% imshow(OtsuMask, []);
% subplot 325;
% imshow(PIMGCropped, []);
% hold on;
% [I, J] = find(OtsuMask);
% plot(J, I, '*');
ExclusionWeight = min(FIMGCropped * 50, 1);
StartMask = OtsuMask & (PIMGCropped > 0.5);
OriginalPSI = double(StartMask);
MeanIntensityInStartMask = mean(ResampledAVWCropped(StartMask));
VarIntensityInStartMask = var(ResampledAVWCropped(StartMask));

P = normpdf(ResampledAVWCropped, MeanIntensityInStartMask, sqrt(VarIntensityInStartMask));
P = P ./ max(P(:));
% subplot 324;
% imshow(P, []);
% subplot 326;

S = OtsuMask & (PIMGCropped > 0.2) & (P > 0.3 | ResampledAVWCropped > MeanIntensityInStartMask);
%StartMask;


%T = S + FIMGCropped;
%imshow(T, []);
% subplot 321;
% hold on;
% [L, CC] = contour(S, [0.5, 0.5]);
%keyboard;
DMask = (P > 0.2 | ResampledAVWCropped > MeanIntensityInStartMask) & ExclusionWeight < 0.5;
% [L, CCT] = contour(DMask, [0.5, 0.5]);
% set(CCT, 'Color', 'r');

S = double(DMask);
S(DMask) = 1;
S(~DMask) = -1;
D = computeDistanceFunction2d(S, [1, 1], [], 2);

OriginalPSI(StartMask) = 1;
OriginalPSI(~StartMask) = -1;
OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);
%ExclusionWeight = log(FIMGCropped + 1);

%FIMGCropped(isinf(FIMGCropped)) = 0;
[PSI, NumIter] = fls_evolution_illusory_contours(-D, OriginalPSI, SIMGCropped, double(ExclusionWeight), [1, 1], 1000, 50, 5, 5, ResampledAVWCropped);
return;

%%
clf;

imshow(ExclusionWeight, []);
hold on;
contour(double(StartMask), [0.5, 0.5]);
%%

%subplot 221;
%hold on;
%[I, J] = find(StartMask);
%plot(J, I, '*');
% hold on;
% imagesc(EIMG, 'AlphaData', AlphaData);
%%
return;
%%
% subplot 131;
% imshow(ResampledAVW, []);
% subplot 132;
% imshow(TemplateAVW, []);
% subplot 133;
% imshow(WarpedIMG, []);

while(1)
	subplot 121;
	imshow(TIMG, []);
	subplot 122;
	imshow(ResampledAVW, []);
	ret = input('press a key');
	subplot 121;
	imshow(ResampledAVW, []);
	subplot 122;
	imshow(TemplateAVW, []);
	ret = input('press a key');
end
