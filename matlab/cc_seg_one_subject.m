function [Overlap] = cc_seg_one_subject(InFile, GroundTruthFile, OutputPNG)

if(exist(InFile, 'file') ~= 2)
	error(['Input file ' InFile ' is not a regular file or does not exist']);
end

[NII, AVW] = load_nii(InFile);
[GroundTruthNII, GroundTruthAVW] = load_nii(GroundTruthFile);

[FornixProbNII, FornixProbAVW] = load_nii('all_fornix_prob.nii.gz');
[ProbNII, ProbAVW] = load_nii('all_cc_prob.nii.gz');
[TemplateNII, TemplateAVW] = load_nii('all_msp_mean.nii.gz');

% [SmoothNII, SmoothAVW] = load_nii('all_msp_smooth_areas.nii.gz');
% 
% SmoothAVW = (SmoothAVW > 0);
% SmoothAVW = imerode(SmoothAVW, ones(3, 3));
% F = gaussian_filter_max_1d(1);
% 
% SmoothAVW = double(SmoothAVW);
% 
% SmoothAVW = imfilter(SmoothAVW, F(:), 'same', 'conv', 'replicate');
% SmoothAVW = imfilter(SmoothAVW, F(:)', 'same', 'conv', 'replicate');
% %SmoothAVW
% %keyboard;
% SmoothAVW = SmoothAVW ./ max(SmoothAVW(:));
% SmoothAVW = SmoothAVW .* 3 + 1;
%[NII, AVW] = load_nii('0001_acpc.nii.gz');

%keyboard;
%[RegNII, RegAVW] = load_nii('0001_acpc_reg.nii.gz');
TemplatePixdims = TemplateNII.hdr.dime.pixdim(3:4);

if ndims(AVW) == 3
	NIIPixdims = NII.hdr.dime.pixdim(3:4);

	MidsagSliceIDX = (size(AVW, 2) + 1) / 2;
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
elseif ndims(AVW) == 2
	AVWMidSag = double(AVW);
	T = NII.hdr.dime.pixdim(2:4);
	NIIPixdims = T(NII.hdr.dime.dim(2:4) > 1);
	clear T;
end

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

R = [163.7248  168.6860  201.4264  125.0233];
%R = [187.3418  195.7478  160.4776   65.7194];
%R = [184.1436  189.6490  165.1613   92.8049];
MatchedFilter = imcrop(TemplateAVW, R);
MatchedFilterProb = imcrop(ProbAVW, R);
MatchedFilterFornixProb = imcrop(FornixProbAVW, R);
%MatchedFilterSmoothAVW = imcrop(SmoothAVW, R);

LeftLimit = max(min(Templatexx), min(AVWMidSagxx));
UpperLimit = max(min(Templateyy), min(AVWMidSagyy));
RightLimit = min(max(Templatexx), max(AVWMidSagxx));
LowerLimit = min(max(Templateyy), max(AVWMidSagyy));

[ResampleX, ResampleY] = meshgrid(LeftLimit:TemplatePixdims(1):RightLimit, UpperLimit:TemplatePixdims(2):LowerLimit);

ResampledAVW = interp2(AVWMidSagX, AVWMidSagY, AVWMidSag, ResampleX, ResampleY);

[THRESH, SEG] = robust_otsu2(ResampledAVW(ResampledAVW > 0), [0.05, 0.98]);
FullSEG = zeros(size(ResampledAVW));
FullSEG(ResampledAVW > 0) = SEG;
FullSEG(ResampledAVW > THRESH(end)) = 3;
FullOtsuSeg = FullSEG;
WMSeg = (FullSEG == 3);

[I, J] = find(WMSeg);

%
%150 pixels * 0.5mm = 7.5cm, 75mm corpus callosum is expected to be 5cm below the top of the skull and in the middle of the image
% we may change this when we start using acpcdetect
% when normcorrx2 runs
% the index cc(I, J) represents the filter's bottom-right point being ResampledAVW(I, J)
TopWMRow = min(I);
% expected centre of the CC is 
ExpectedCentreRow = round(TopWMRow + 75 ./ TemplatePixdims(2));
ExpectedCentreCol = round(size(ResampledAVW, 2) / 2);

ExpectedLowerRightI = round(ExpectedCentreRow + size(MatchedFilter, 1) / 2);
ExpectedLowerRightJ = round(ExpectedCentreCol + size(MatchedFilter, 2) / 2);

F = gaussian_filter_max_1d(2);
WMSegSmoothed = imfilter(double(WMSeg), F(:), 'same', 'conv', 'replicate');
WMSegSmoothed = imfilter(WMSegSmoothed, F(:)', 'same', 'conv', 'replicate');

cc = normxcorr2(MatchedFilter, ResampledAVW);
cc2 = normxcorr2(MatchedFilterProb, WMSegSmoothed);

[cc2X, cc2Y] = meshgrid((((1:size(cc2, 2)) - ExpectedLowerRightJ)) .* TemplatePixdims(1), (((1:size(cc2, 1)) - ExpectedLowerRightI)) .* TemplatePixdims(2));

R = sqrt(cc2X .* cc2X + cc2Y .* cc2Y);
RWeighting = 0.5 * tanh(-(R - 50) ./ 10) + 0.5; 
cc2R = cc2 .* RWeighting;
clear SEG;

% subplot 121;
% imshow(ResampledAVW, []);
% 
% subplot 122;
% imshow(RWeighting, []);
% keyboard;
%L = bwlabel(C == 3);
%keyboard;
%%
%MatchedFilterProb
% subplot 221;
% imshow(cc, []);
% subplot 222;
% imshow(cc2, []);
% subplot 223;
% imshow(CSmoothed, []);
% %%
%keyboard;

[~, SortedCCIDX] = sort(cc(:), 1, 'descend');

[SortedCCID, SortedCCJD] = ind2sub(size(cc), SortedCCIDX);

yoffset = round(SortedCCID - size(MatchedFilter, 1));
xoffset = round(SortedCCJD - size(MatchedFilter, 2));

FirstOffsetValid = find(yoffset >= 1 & yoffset + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2) & xoffset >= 1 & xoffset + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1), 1, 'first');

total_offset = [yoffset(FirstOffsetValid), xoffset(FirstOffsetValid)];
corr_offset = total_offset;
xpeak = xoffset(FirstOffsetValid);
ypeak = yoffset(FirstOffsetValid);

%keyboard;
[~, SortedCCIDX] = sort(cc2(:), 1, 'descend');

[SortedCCID, SortedCCJD] = ind2sub(size(cc2), SortedCCIDX);

yoffset2 = round(SortedCCID - size(MatchedFilter, 1));
xoffset2 = round(SortedCCJD - size(MatchedFilter, 2));

FirstOffsetValid2 = find(yoffset2 >= 1 & yoffset2 + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2) & xoffset2 >= 1 & xoffset2 + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1), 1, 'first');

total_offset2 = [yoffset2(FirstOffsetValid2), xoffset2(FirstOffsetValid2)];
corr_offset2 = total_offset2;
xpeak2 = xoffset2(FirstOffsetValid2);
ypeak2 = yoffset2(FirstOffsetValid2);

[~, SortedCCIDX] = sort(cc2R(:), 1, 'descend');

[SortedCCID, SortedCCJD] = ind2sub(size(cc2R), SortedCCIDX);

yoffset2R = round(SortedCCID - size(MatchedFilter, 1));
xoffset2R = round(SortedCCJD - size(MatchedFilter, 2));

FirstOffsetValid2R = find(yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2) & xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1), 1, 'first');

total_offset2R = [yoffset2R(FirstOffsetValid2R), xoffset2R(FirstOffsetValid2R)];
corr_offset2R = total_offset2R;
xpeak2R = xoffset2R(FirstOffsetValid2R);
ypeak2R = yoffset2R(FirstOffsetValid2R);

real_total_offset = total_offset2R;
real_corr_offset = corr_offset2R;
realxpeak = xpeak2R;
realypeak = ypeak2R;
%keyboard;
%[max_cc, imax] = max(abs(cc(:)));

%[ypeak, xpeak] = ind2sub(size(cc),imax(1));
%corr_offset_orig = [ (ypeak-size(MatchedFilter,1)) (xpeak-size(MatchedFilter,2)) ];
%keyboard;

%corr_offset = round(corr_offset + R([2 1]));

%rect_offset = [size(ResampledAVW, 1) - size(MatchedFilter, 1), size(ResampledAVW, 2) - size(MatchedFilter, 2)];

%total_offset = round(corr_offset);
%E = edge(MatchedFilter, 'canny', [0.1 0.4]);
%E(:) = 1;
%EIMG = repmat(reshape([0, 0, 1], [1, 1, 3]), [size(ResampledAVW), 1]);
%AlphaData = zeros(size(ResampledAVW));
%AlphaData = E;
%F = ResampledAVW;
%AlphaData(total_offset(1):total_offset(1) + size(MatchedFilter, 1) - 1, total_offset(2):total_offset(2) + size(MatchedFilter, 2) - 1) = E;

%keyboard;
% do histogram matching
InitialResampledAVW = ResampledAVW(real_total_offset(1):real_total_offset(1) + size(MatchedFilter, 1) - 1, real_total_offset(2):real_total_offset(2) + size(MatchedFilter, 2) - 1);

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
%disp('Running LK');

%I = (ResampledAVW - min(ResampledAVW(:))) ./ (max(ResampledAVW(:)) - min(ResampledAVW(:)));
%J = (MatchedFilter - min(MatchedFilter(:))) ./ (max(MatchedFilter(:)) - min(MatchedFilter(:)));
Parameters = lk_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, [0 0 0 0 corr_offset2R(2) corr_offset2R(1)], 100);
ParametersW = lk_weighted_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, MatchedFilter, [0 0 0 0 corr_offset2R(2) corr_offset2R(1)], 100);

T = MatchedFilter;
T = T ./ max(T(:));
ParametersSeg = lk_run_affine_inv_comp(WMSegSmoothed, T, [0 0 0 0 corr_offset2R(2) corr_offset2R(1)], 100);
ParametersSegW = lk_weighted_run_affine_inv_comp(WMSegSmoothed, T, MatchedFilter, [0 0 0 0 corr_offset2R(2) corr_offset2R(1)], 100);

[TX, TY, InterpX, InterpY] = coords_template_lk_img(Parameters, ResampledAVW, MatchedFilter);
[TXW, TYW, InterpXW, InterpYW] = coords_template_lk_img(ParametersW, ResampledAVW, MatchedFilter);
[TXSeg, TYSeg, InterpXSeg, InterpYSeg] = coords_template_lk_img(ParametersSeg, ResampledAVW, MatchedFilter);
[TXSegW, TYSegW, InterpXSegW, InterpYSegW] = coords_template_lk_img(ParametersSeg, ResampledAVW, MatchedFilter);

% [TemplateXISO, TemplateYISO] = meshgrid(1:size(MatchedFilter, 2), 1:size(MatchedFilter, 1));
% %TemplateXISO = TemplateXISO - corr_offset(2);
% %TemplateYISO = TemplateYISO - corr_offset(1);
% XY = cat(1, TemplateXISO(:)', TemplateYISO(:)', ones(1, numel(MatchedFilter)));
% %%
% TransformationMatrix = [1 + Parameters(1), Parameters(3), Parameters(5); ...
%     Parameters(2), 1 + Parameters(4), Parameters(6); ...
%     0, 0, 1];
% 
% TXY = TransformationMatrix * XY;
% 
% TXY = TXY(1:2, :);
% 
% TX = reshape(TXY(1, :), size(MatchedFilter));
% TY = reshape(TXY(2, :), size(MatchedFilter));
% clear TXY TransformationMatrix;

%% weighted version
% TransformationMatrixW = [1 + ParametersW(1), ParametersW(3), ParametersW(5); ...
%     ParametersW(2), 1 + ParametersW(4), ParametersW(6); ...
%     0, 0, 1];
% 
% TXYW = TransformationMatrixW * XY;
% 
% TXYW = TXYW(1:2, :);
% 
% TXW = reshape(TXYW(1, :), size(MatchedFilter));
% TYW = reshape(TXYW(2, :), size(MatchedFilter));
% clear TXYW TransformationMatrixW;
% 
% %% seg version
% TransformationMatrixSegW = [1 + ParametersSegW(1), ParametersSegW(3), ParametersSegW(5); ...
%     ParametersSegW(2), 1 + ParametersSegW(4), ParametersSegW(6); ...
%     0, 0, 1];
% 
% TXYSegW = TransformationMatrixSegW * XY;
% 
% TXYSegW = TXYSegW(1:2, :);
% 
% TXSegW = reshape(TXYSegW(1, :), size(MatchedFilter));
% TYSegW = reshape(TXYSegW(2, :), size(MatchedFilter));
% clear TXYSegW TransformationMatrixSegW;
% 
% %% seg version
% TransformationMatrixSeg = [1 + ParametersSeg(1), ParametersSeg(3), ParametersSeg(5); ...
%     ParametersSeg(2), 1 + ParametersSeg(4), ParametersSeg(6); ...
%     0, 0, 1];
% 
% TXYSeg = TransformationMatrixSeg * XY;
% 
% TXYSeg = TXYSeg(1:2, :);
% 
% TXSeg = reshape(TXYSeg(1, :), size(MatchedFilter));
% TYSeg = reshape(TXYSeg(2, :), size(MatchedFilter));
% clear TXYSeg TransformationMatrixSeg;
% 
% % find the change in vector along the X direction
% XChangeVector = [TX(1, 2) - TX(1, 1), TY(1, 2) - TY(1, 1)];
% YChangeVector = [TX(2, 1) - TX(1, 1), TY(2, 1) - TY(1, 1)];
% 
% UnitXChangeVector = XChangeVector ./ norm(XChangeVector);
% UnitYChangeVector = YChangeVector ./ norm(YChangeVector);	
% 
% % find the change in vector along the X direction
% % scalar resolute in the direction of XChangeVector
% 
% [ResampledAVWX, ResampledAVWY] = meshgrid(1:size(ResampledAVW, 2), 1:size(ResampledAVW, 1));
% 
% XDot = (ResampledAVWX - TX(1, 1)) .* UnitXChangeVector(1) + (ResampledAVWY - TY(1, 1)) .* UnitXChangeVector(2);
% YDot = (ResampledAVWX - TX(1, 1)) .* UnitYChangeVector(1) + (ResampledAVWY - TY(1, 1)) .* UnitYChangeVector(2);
% 
% X = XDot ./ norm(XChangeVector);
% Y = YDot ./ norm(YChangeVector);
% % interpolate 
%WarpedIMG = interp2(TemplateAVW, TX, TY, 'bilinear');

TIMG = interp2(MatchedFilter, InterpX, InterpY, 'linear', 0);
PIMG = interp2(MatchedFilterProb, InterpX, InterpY, 'linear', 0);
FIMG = interp2(MatchedFilterFornixProb, InterpX, InterpY, 'linear', 0);
PIMGW = interp2(MatchedFilterProb, InterpXW, InterpYW, 'linear', 0);
PIMGSeg = interp2(MatchedFilterProb, InterpXSeg, InterpYSeg, 'linear', 0);
TIMGSegW = interp2(MatchedFilter, InterpXSegW, InterpYSegW, 'linear', 0);
PIMGSegW = interp2(MatchedFilterProb, InterpXSegW, InterpYSegW, 'linear', 0);

CCSeg = (GroundTruthAVW > 0);
Overlap.Original = sum(PIMG(:) .* CCSeg(:)) ./ sum(PIMG(:));
Overlap.OriginalW = sum(PIMGW(:) .* CCSeg(:)) ./ sum(PIMGW(:));
Overlap.Seg = sum(PIMGSeg(:) .* CCSeg(:)) ./ sum(PIMGSeg(:));
Overlap.SegW = sum(PIMGSegW(:) .* CCSeg(:)) ./ sum(PIMGSegW(:));
return;
%disp('
% %%
% imshow(TIMGSegW, [])
% hold on;
% plot(TXSegW, TYSegW);
%%
%keyboard;
%SIMG = interp2(MatchedFilterSmoothAVW, X, Y, 'linear', 0);
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
if(length(R) > 1)
	disp('Bad image, returning');
	return;
end
ResampledAVWCropped = imcrop(ResampledAVW, R.BoundingBox);
PIMGCropped = imcrop(PIMG, R.BoundingBox);
FIMGCropped = imcrop(FIMG, R.BoundingBox);
%SIMGCropped = imcrop(SIMG, R.BoundingBox);

[THRESH, SEG] = robust_otsu2(ResampledAVWCropped, [0.05, 0.98]);
SEG(ResampledAVWCropped > THRESH(end)) = 3;
%keyboard;
ResampledAVWCroppedForOtsu = (ResampledAVWCropped - min(ResampledAVWCropped(:))) ./ (max(ResampledAVWCropped(:)) - min(ResampledAVWCropped(:)));
ResampledAVWCroppedForOtsu = round(ResampledAVWCroppedForOtsu * 255);
%%

%OriginalOtsuSeg = reshape(OriginalOtsuSeg, size(ResampledAVWForOtsu));
%keyboard;

T = otsu2(ResampledAVWCroppedForOtsu);
OtsuMask = (ResampledAVWCroppedForOtsu >= T(end));
OriginalOtsuMask = OtsuMask;

% do a second LK on the probability map
%%


%%
TemplateOverlapOtsu = sum(OtsuMask(:) .* PIMGCropped(:)) ./ sum(PIMGCropped(:));
[OriginalOtsuSegCounts, OriginalOtsuSeg] = histc(ResampledAVWCroppedForOtsu(:), [0, T, 256]);
OriginalOtsuSeg = reshape(OriginalOtsuSeg, size(ResampledAVWCroppedForOtsu));
%keyboard;
OtsuMask = imopen(OtsuMask, strel('diamond', 2));
%%
% subplot 121;
% imshow(S, []);
% subplot 122;
% imshow(OtsuMask, []);
% 
% %%
% keyboard;

FornixMaskIdx = find((FIMGCropped * 50) > 0.5);

L = bwlabel(OtsuMask);
R = regionprops(L, 'Area', 'PixelIdxList');

TemplateOverlap = zeros(length(R), 1);
% keep all regions that have high overlap with template or big area
for z = 1:length(R)
	IDX = setdiff(R(z).PixelIdxList, FornixMaskIdx);
	TemplateOverlap(z) = sum(PIMGCropped(IDX)) ./ length(IDX);
	clear IDX;
end
clear FornixMaskIdx;
%keyboard;
%[~, I] = max([R.Area]);

OtsuMask = ismember(L, find(TemplateOverlap > 0.4));
clear Junk;

ExclusionWeight = min(FIMGCropped * 50, 1);

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

StartMask = OtsuMask & (PIMGCropped > 0.5);
OriginalPSI = double(StartMask);
MeanIntensityInStartMask = mean(ResampledAVWCropped(StartMask));
VarIntensityInStartMask = var(ResampledAVWCropped(StartMask));

GaussianProb = normpdf(ResampledAVWCropped, MeanIntensityInStartMask, sqrt(VarIntensityInStartMask));
GaussianProb = GaussianProb ./ max(GaussianProb(:));
GaussianProb(ResampledAVWCropped > MeanIntensityInStartMask) = 1;

SIGMA = 0.5;

[GaussianFilter, GaussianFilterDeriv] = gaussian_filter_max_1d(SIGMA);

aSmooth=imfilter(GaussianProb, GaussianFilter, 'conv', 'replicate', 'same');   % run the filter across rows
aSmooth=imfilter(aSmooth, GaussianFilter','conv','replicate', 'same'); % and then across columns

%apply directional derivatives
ax = imfilter(aSmooth, GaussianFilterDeriv(:)', 'conv','replicate');
ay = imfilter(aSmooth, GaussianFilterDeriv(:), 'conv','replicate');

GaussianGradMAG = sqrt(ax .* ax + ay .* ay);
GaussianGradMAG = GaussianGradMAG ./ max(GaussianGradMAG(:));
ExclusionWeight = ExclusionWeight + GaussianGradMAG + exp(-GaussianProb);
% subplot 324;
% imshow(P, []);
% subplot 326;

S = OtsuMask & (PIMGCropped > 0.2) & (GaussianProb > 0.3 | ResampledAVWCropped > MeanIntensityInStartMask);
%StartMask;

%T = S + FIMGCropped;
%imshow(T, []);
% subplot 321;
% hold on;
% [L, CC] = contour(S, [0.5, 0.5]);
%keyboard;

DMask = (GaussianProb > 0.2 | ResampledAVWCropped > MeanIntensityInStartMask) & ExclusionWeight < 0.5;
DMask = imdilate(DMask, ones(5, 5));
% [L, CCT] = contour(DMask, [0.5, 0.5]);
% set(CCT, 'Color', 'r');

S = double(DMask);
S(DMask) = 1;
S(~DMask) = -1;
D = computeDistanceFunction2d(S, [1, 1], [], 2);

OriginalPSI(StartMask) = 1;
OriginalPSI(~StartMask) = -1;
%OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);
OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);
%OriginalPSI = fls_solve_eikonal_equation(OriginalPSI, [1, 1], 'ENO1', true(size(OriginalPSI)), Inf, 2);
%%
%subplot 121; imshow(OriginalPSILSMLIB, []);
%subplot 122; imshow(OriginalPSI, []);
%%
%keyboard;
%CurvatureDSecondOrder = ls_curvature(D, H, true(size(D)), 'gaussian');
%CurvatureDFourthOrder = ls_curvature_fourth_order(D, H, true(size(D)), 'gaussian');
%ExclusionWeight = log(FIMGCropped + 1);
% 
% clf;
% subplot 121;
% imshow(-D > 0, []);
% hold on;
% [C, CC] = contour(OriginalPSI, [0, 0]);
% set(CC, 'Color', 'r');
% subplot 122;
% imshow(OriginalPSI > 0, []);
% hold on;
% [C, CC] = contour(OriginalPSI, [0, 0]);
% set(CC, 'Color', 'r');
% 
%keyboard;
%if(all(D(:) > 0) || all(D(:) < 0) || any(isnan(GaussianProb(:))))
ContourLevel = 0.15;
	clf;
	subplot 221;
	imshow(ResampledAVW, []);
	hold on;
	SZ = size(MatchedFilter);

	T = line([xpeak, xpeak + SZ(2)], [ypeak, ypeak]); set(T, 'Color', 'b');
	T = line([xpeak, xpeak + SZ(2)], [ypeak + SZ(1), ypeak + SZ(1)]); set(T, 'Color', 'b');
	T = line([xpeak, xpeak], [ypeak, ypeak + SZ(1)]); set(T, 'Color', 'b');
	T = line([xpeak + SZ(2), xpeak + SZ(2)], [ypeak, ypeak + SZ(1)]); set(T, 'Color', 'b');

	T = line([xpeak2, xpeak2 + SZ(2)], [ypeak2, ypeak2]); set(T, 'Color', 'g');
	T = line([xpeak2, xpeak2 + SZ(2)], [ypeak2 + SZ(1), ypeak2 + SZ(1)]); set(T, 'Color', 'g');
	T = line([xpeak2, xpeak2], [ypeak2, ypeak2 + SZ(1)]); set(T, 'Color', 'g');
	T = line([xpeak2 + SZ(2), xpeak2 + SZ(2)], [ypeak2, ypeak2 + SZ(1)]); set(T, 'Color', 'g');

	T = line([xpeak2R, xpeak2R + SZ(2)], [ypeak2R, ypeak2R]); set(T, 'Color', 'm');
	T = line([xpeak2R, xpeak2R + SZ(2)], [ypeak2R + SZ(1), ypeak2R + SZ(1)]); set(T, 'Color', 'm');
	T = line([xpeak2R, xpeak2R], [ypeak2R, ypeak2R + SZ(1)]); set(T, 'Color', 'm');
	T = line([xpeak2R + SZ(2), xpeak2R + SZ(2)], [ypeak2R, ypeak2R + SZ(1)]); set(T, 'Color', 'm');

	%plot(TX(:, 1), TY(:, 1), 'r*');
	%plot(TX(:, end), TY(:, end), 'r*');
	%plot(TX(1, :), TY(1, :), 'r*');
	%plot(TX(end, :), TY(end, :), 'r*');
	
	%plot(TX, TY, '*');
	title('Initial XCORR alignment, blue entire image, green seg, magenta seg with weight');
	subplot 222;
	imshow(ResampledAVW, []);
	hold on;
	SZ = size(MatchedFilter);
	
	T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
	T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
	T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
	T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
	
	T = line([TXW(1, 1) TXW(end, 1)], [TYW(1, 1), TYW(end, 1)]); set(T, 'Color', 'y');
	T = line([TXW(1, end) TXW(end, end)], [TYW(1, end), TYW(end, end)]); set(T, 'Color', 'y');
	T = line([TXW(1, 1) TXW(1, end)], [TYW(1, 1), TYW(1, end)]); set(T, 'Color', 'y');
	T = line([TXW(end, 1) TXW(end, end)], [TYW(end, 1), TYW(end, end)]); set(T, 'Color', 'y');
	[~, CC] = contour(PIMG, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
	[~, CC] = contour(PIMGW, [ContourLevel, ContourLevel]); set(CC, 'Color', 'y');
	
	title('Initial LK alignment, red uniform, yellow weighted');
	
	subplot 223;
	imshow(ResampledAVW, []);
	hold on;
	SZ = size(MatchedFilter);
	
	T = line([TXSeg(1, 1) TXSeg(end, 1)], [TYSeg(1, 1), TYSeg(end, 1)]); set(T, 'Color', 'r');
	T = line([TXSeg(1, end) TXSeg(end, end)], [TYSeg(1, end), TYSeg(end, end)]); set(T, 'Color', 'r');
	T = line([TXSeg(1, 1) TXSeg(1, end)], [TYSeg(1, 1), TYSeg(1, end)]); set(T, 'Color', 'r');
	T = line([TXSeg(end, 1) TXSeg(end, end)], [TYSeg(end, 1), TYSeg(end, end)]); set(T, 'Color', 'r');
	
	T = line([TXSegW(1, 1) TXSegW(end, 1)], [TYSegW(1, 1), TYSegW(end, 1)]); set(T, 'Color', 'y');
	T = line([TXSegW(1, end) TXSegW(end, end)], [TYSegW(1, end), TYSegW(end, end)]); set(T, 'Color', 'y');
	T = line([TXSegW(1, 1) TXSegW(1, end)], [TYSegW(1, 1), TYSegW(1, end)]); set(T, 'Color', 'y');
	T = line([TXSegW(end, 1) TXSegW(end, end)], [TYSegW(end, 1), TYSegW(end, end)]); set(T, 'Color', 'y');

	[~, CC] = contour(PIMGSeg, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
	[~, CC] = contour(PIMGSegW, [ContourLevel, ContourLevel]); set(CC, 'Color', 'y');
	title('Initial LK alignment to seg, red uniform, yellow weighted');
	
% 	subplot 223;
% 	%RWeightingIMG = padarray(RWeighting, size(MatchedFilter), 0, 'pre');
% 	imshow(cc2R, []);
% 	title('Correlation coefficient output');
	subplot 224;
	%RWeightingIMG = padarray(RWeighting, size(MatchedFilter), 0, 'pre');
	imshow(PIMG, []);
	hold on;
	T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
	T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
	T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
	T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
	
	[~, CC] = contour(PIMG, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
	
	%imshow(WMSegSmoothed, []);
	%title('WMSeg Smoothed');
	title('P');
	FigPos = fullscreen_fig_pos;
	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
	exportfig(gcf, [OutputPNG '-cc2'], 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

	disp('Bad image, returning');
	return;
%else
%	return;
%end
return;
%keyboard;
%FIMGCropped(isinf(FIMGCropped)) = 0;
[NewPSI, NumIter] = fls_evolution_illusory_contours(-D, OriginalPSI, SIMGCropped, double(ExclusionWeight), [1, 1], 1000, 200, 2, 5, ResampledAVWCropped);
clf;
SR = 3;
SC = 3;
subplot(SR, SC, 1);
imshow(ResampledAVW, []);
hold on;
SZ = size(MatchedFilter);

T = line([xpeak, xpeak + SZ(2)], [ypeak, ypeak]); set(T, 'Color', 'b');
T = line([xpeak, xpeak + SZ(2)], [ypeak + SZ(1), ypeak + SZ(1)]); set(T, 'Color', 'b');
T = line([xpeak, xpeak], [ypeak, ypeak + SZ(1)]); set(T, 'Color', 'b');
T = line([xpeak + SZ(2), xpeak + SZ(2)], [ypeak, ypeak + SZ(1)]); set(T, 'Color', 'b');
T = line([xpeak2, xpeak2 + SZ(2)], [ypeak2, ypeak2]); set(T, 'Color', 'g');
T = line([xpeak2, xpeak2 + SZ(2)], [ypeak2 + SZ(1), ypeak2 + SZ(1)]); set(T, 'Color', 'g');
T = line([xpeak2, xpeak2], [ypeak2, ypeak2 + SZ(1)]); set(T, 'Color', 'g');
T = line([xpeak2 + SZ(2), xpeak2 + SZ(2)], [ypeak2, ypeak2 + SZ(1)]); set(T, 'Color', 'g');

%plot(TX(:, 1), TYSeg(:, 1), 'r*');
%plot(TX(:, end), TY(:, end), 'r*');
%plot(TX(1, :), TY(1, :), 'r*');
%plot(TX(end, :), TY(end, :), 'r*');
T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
%plot(TX, TY, '*');
title('Initial XCORR and LK alignment');

subplot(SR, SC, 2);
imshow(ResampledAVWCropped, []);
hold on;
[C, CC] = contour(NewPSI, [0, 0]);
set(CC, 'Color', 'r');SortedIMGTruncForOtsu
[C, CC] = contour(OriginalPSI, [0, 0]);
set(CC, 'Color', 'b');
title('Original Image');

subplot(SR, SC, 3);
imshow(GaussianProb, []);
hold on;
[C, CC] = contour(NewPSI, [0, 0]);
set(CC, 'Color', 'r');
[C, CC] = contour(OriginalPSI, [0, 0]);
set(CC, 'Color', 'b');
title('Gaussian probability of pixels');

subplot(SR, SC, 4);
imshow(PIMGCropped, []);
hold on;
[C, CC] = contour(NewPSI, [0, 0]);
set(CC, 'Color', 'r');
[C, CC] = contour(OriginalPSI, [0, 0]);
set(CC, 'Color', 'b');
title('Probability of CC from template');

subplot(SR, SC, 5);
imshow(PIMGCropped, []);
hold on;
[C, CC] = contour(OriginalOtsuMask, [0.5, 0.5]);
set(CC, 'Color', 'r');
title(['Probability of CC from template with original otsu borders, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']);

subplot(SR, SC, 6);
imshow(OriginalOtsuSeg, []);
hold on;
[C, CC] = contour(NewPSI, [0, 0]);
set(CC, 'Color', 'r');
[C, CC] = contour(OriginalPSI, [0, 0]);
set(CC, 'Color', 'b');
title(['OtsuMask ' num2str(OriginalOtsuSegCounts(1:3)', '%d ')]);

subplot(SR, SC, 7);
imshow(ExclusionWeight, []);
hold on;
[C, CC] = contour(NewPSI, [0, 0]);
set(CC, 'Color', 'r');
[C, CC] = contour(OriginalPSI, [0, 0]);
set(CC, 'Color', 'b');
title('Exclusion mask');

FigPos = fullscreen_fig_pos;
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

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

function [TX, TY, InterpX, InterpY] = coords_template_lk_img(Parameters, IMG, Template)

% returns the coordinates of thetemplate transformed by the LK affine transformation matrix Parameters in the space of the target image IMG

[TemplateX, TemplateY] = meshgrid(1:size(Template, 2), 1:size(Template, 1));

XY = cat(1, TemplateX(:)', TemplateY(:)', ones(1, numel(Template)));

TransformationMatrix = [1 + Parameters(1), Parameters(3), Parameters(5); ...
    Parameters(2), 1 + Parameters(4), Parameters(6); ...
    0, 0, 1];

TXY = TransformationMatrix * XY;

TXY = TXY(1:2, :);

TX = reshape(TXY(1, :), size(Template));
TY = reshape(TXY(2, :), size(Template));
clear TXY;

% dot product method doesnt seem to work
% % find the change in vector along the X direction
% XChangeVector = [TX(1, 2) - TX(1, 1), TY(1, 2) - TY(1, 1)];
% % find the change in vector along the Y direction
% YChangeVector = [TX(2, 1) - TX(1, 1), TY(2, 1) - TY(1, 1)];
% 
% UnitXChangeVector = XChangeVector ./ norm(XChangeVector);
% UnitYChangeVector = YChangeVector ./ norm(YChangeVector);	
% 
% % find the change in vector along the X direction
% % scalar resolute in the direction of XChangeVector
% 
% [IMGX, IMGY] = meshgrid(1:size(IMG, 2), 1:size(IMG, 1));
% 
% XDot = (IMGX - TX(1, 1)) .* UnitXChangeVector(1) + (IMGY - TY(1, 1)) .* UnitXChangeVector(2);
% YDot = (IMGX - TX(1, 1)) .* UnitYChangeVector(1) + (IMGY - TY(1, 1)) .* UnitYChangeVector(2);
% 
% InterpX = XDot ./ norm(XChangeVector);
% InterpY = YDot ./ norm(YChangeVector);

[IMGX, IMGY] = meshgrid(1:size(IMG, 2), 1:size(IMG, 1));

ITXY = TransformationMatrix \ [IMGX(:)'; IMGY(:)'; ones(1, numel(IMGX))];
InterpX = reshape(ITXY(1, :), size(IMGX));
InterpY = reshape(ITXY(2, :), size(IMGX));
% keyboard;
% 
% %%
% clf;
% subplot 221;
% imshow(IMG, []);
% hold on;
% plot(TX, TY, '*');
% IMGI = interp2(Template, InterpX, InterpY, 'linear', 0);
% subplot 222;
% imshow(IMGI, []);
% T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
% T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
% T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
% T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
% IMGI = interp2(Template, ITX, ITY, 'linear', 0);
% subplot 223;
% imshow(IMGI, []);
% T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
% T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
% T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
% T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
