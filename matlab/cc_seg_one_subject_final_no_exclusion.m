function [UsedConvexHull] = cc_seg_one_subject_final_no_exclusion(InFile, GroundTruthFile, OutputDir, OutputPrefix)

if(exist(InFile, 'file') ~= 2)
	error(['Input file ' InFile ' is not a regular file or does not exist']);
end

[~, ~, ~] = mkdir(OutputDir);

OutputMAT = fullfile(OutputDir, [OutputPrefix '.mat']);
OutputPNG = fullfile(OutputDir, [OutputPrefix '.png']);

[NII, AVW] = load_nii(InFile);
[GroundTruthNII, GroundTruthAVW] = load_nii(GroundTruthFile);

[FornixProbNII, FornixProbAVW] = load_nii('all_fornix_prob.nii.gz');
[ProbNII, ProbAVW] = load_nii('all_cc_prob.nii.gz');
[TemplateNII, TemplateAVW] = load_nii('all_msp_mean.nii.gz');
[SmoothNII, SmoothAVW] = load_nii('all_msp_smooth_areas.nii.gz');

%DilationANGLES = load('all_cc_dilation_angle');
% 
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

%[TemplateX, TemplateY] = meshgrid(Templatexx, Templateyy);

R = [163.7248  168.6860  201.4264  125.0233];
%R = [187.3418  195.7478  160.4776   65.7194];
%R = [184.1436  189.6490  165.1613   92.8049];
MatchedFilter = imcrop(TemplateAVW, R);
MatchedFilterProb = imcrop(ProbAVW, R);
MatchedFilterFornixProb = imcrop(FornixProbAVW, R);
MatchedFilterSmoothAVW = imcrop(SmoothAVW, R);
%MatchedFilterAngles = imcrop(DilationANGLES.ANGLE, R);

LeftLimit = max(min(Templatexx), min(AVWMidSagxx));
UpperLimit = max(min(Templateyy), min(AVWMidSagyy));
RightLimit = min(max(Templatexx), max(AVWMidSagxx));
LowerLimit = min(max(Templateyy), max(AVWMidSagyy));

[ResampleX, ResampleY] = meshgrid(LeftLimit:TemplatePixdims(1):RightLimit, UpperLimit:TemplatePixdims(2):LowerLimit);

ResampledAVW = interp2(AVWMidSagX, AVWMidSagY, AVWMidSag, ResampleX, ResampleY);
ResampledGroundAVW = interp2(AVWMidSagX, AVWMidSagY, GroundTruthAVW, ResampleX, ResampleY, 'nearest');

[THRESH, SEG] = robust_otsu2(ResampledAVW(ResampledAVW > 0), [0.05, 0.98]);
FullSEG = zeros(size(ResampledAVW));
FullSEG(ResampledAVW > 0) = SEG;
FullSEG(ResampledAVW > THRESH(end)) = 3;
%FullOtsuSeg = FullSEG;
WMSeg = (FullSEG == 3);

[I, J] = find(WMSeg);
clear J;
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

cc2 = normxcorr2(MatchedFilterProb, WMSegSmoothed);

[cc2X, cc2Y] = meshgrid((((1:size(cc2, 2)) - ExpectedLowerRightJ)) .* TemplatePixdims(1), (((1:size(cc2, 1)) - ExpectedLowerRightI)) .* TemplatePixdims(2));

R = sqrt(cc2X .* cc2X + cc2Y .* cc2Y);
RWeighting = 0.5 * tanh(-(R - 50) ./ 10) + 0.5; 
cc2R = cc2 .* RWeighting;
clear SEG;

cc2RRegMax = imregionalmax(cc2R, 4) & (cc2R > 0);
[cc2RRegMaxI, cc2RRegMaxJ] = find(cc2RRegMax);
Rofcc2RRegMax = R(cc2RRegMax);
[~, I] = min(Rofcc2RRegMax);

cc2RRegMaxIMaxI = cc2RRegMaxI(I);
cc2RRegMaxIMaxJ = cc2RRegMaxJ(I);
yoffset2RReg = round(cc2RRegMaxI(I) - size(MatchedFilter, 1));
xoffset2RReg = round(cc2RRegMaxJ(I) - size(MatchedFilter, 2));

[~, SortedCCIDX] = sort(cc2R(:), 1, 'descend');
[SortedCCID, SortedCCJD] = ind2sub(size(cc2R), SortedCCIDX);

yoffset2R = round(SortedCCID - size(MatchedFilter, 1));
xoffset2R = round(SortedCCJD - size(MatchedFilter, 2));

FirstOffsetValid2R = find(yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2) & xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1), 1, 'first');

total_offset2R = [yoffset2R(FirstOffsetValid2R), xoffset2R(FirstOffsetValid2R)];
corr_offset2R = total_offset2R;
xpeak2R = xoffset2R(FirstOffsetValid2R);
ypeak2R = yoffset2R(FirstOffsetValid2R);

real_total_offset = [yoffset2RReg, xoffset2RReg];
realxpeak = xoffset2RReg;
realypeak = yoffset2RReg;
%keyboard;
%real_total_offset = total_offset2R;
%real_corr_offset = corr_offset2R;
%realxpeak = xpeak2R;
%realypeak = ypeak2R;
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

MatchedFilterRemapped = immatchhist(MatchedFilter, InitialResampledAVW);

% InitialResampledAVWBins = linspace(min(InitialResampledAVW(:)), max(InitialResampledAVW(:)), 998);
% S = InitialResampledAVWBins(2) - InitialResampledAVWBins(1);
% InitialResampledAVWBins = [InitialResampledAVWBins(1) - S, InitialResampledAVWBins, InitialResampledAVWBins(end) + S];
% clear S;
% 
% % get the CDFs of the histogram of the target image
% [InitialResampledAVWHist, InitialResampledAVWBins] = hist(InitialResampledAVW(:), InitialResampledAVWBins);
% InitialResampledAVWHist = InitialResampledAVWHist ./ sum(InitialResampledAVWHist);
% InitialResampledAVWHistCDF = cumsum(InitialResampledAVWHist);
% % get the CDFs of the histogram of the matched filter
% 
% MatchedFilterBins = linspace(min(MatchedFilter(:)), max(MatchedFilter(:)), 998);
% S = MatchedFilterBins(2) - MatchedFilterBins(1);
% MatchedFilterBins = [MatchedFilterBins(1) - S, MatchedFilterBins, MatchedFilterBins(end) + S];
% clear S;
% [MatchedFilterHist, MatchedFilterBins] = hist(MatchedFilter(:), MatchedFilterBins);
% MatchedFilterHist = MatchedFilterHist ./ sum(MatchedFilterHist);
% MatchedFilterHistCDF = cumsum(MatchedFilterHist);
% 
% % use unique to create a CDF with no repeated elements
% [UMatchedFilterHistCDF, I, J] = unique(MatchedFilterHistCDF);
% UMatchedFilterHistCDFBins = MatchedFilterBins(I);
% 
% [UInitialResampledAVWHistCDF, I, J] = unique(InitialResampledAVWHistCDF);
% UInitialResampledAVWHistCDFBins = InitialResampledAVWBins(I);
% clear I J;
% 
% CDFMatchedFilterPixels = interp1(UMatchedFilterHistCDFBins, UMatchedFilterHistCDF, MatchedFilter(:));
% MatchedFilterRemapped = interp1(UInitialResampledAVWHistCDF, UInitialResampledAVWHistCDFBins, CDFMatchedFilterPixels);
% 
% MatchedFilterRemapped = reshape(MatchedFilterRemapped, size(MatchedFilter));
%keyboard;
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
%Parameters = lk_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, [0 0 0 0 corr_offset2R(2) corr_offset2R(1)], 100);
Parameters = lk_weighted_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, MatchedFilter, [0 0 0 0 real_total_offset(2) real_total_offset(1)], 100);

[TX, TY, InterpX, InterpY] = coords_template_lk_img(Parameters, ResampledAVW, MatchedFilter);

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
SIMG = interp2(MatchedFilterSmoothAVW, InterpX, InterpY, 'linear', 0);
%AIMG = interp2(MatchedFilterAngles, InterpX, InterpY, 'nearest');
%[~, L] = bwdist(~isnan(AIMG));
%AIMG = AIMG(L);
%clear L;
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
ResampledGroundCropped = imcrop(ResampledGroundAVW, R.BoundingBox);
PIMGCropped = imcrop(PIMG, R.BoundingBox);
FIMGCropped = imcrop(FIMG, R.BoundingBox);
SIMGCropped = imcrop(SIMG, R.BoundingBox);
%AIMGCropped = imcrop(AIMG, R.BoundingBox);

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
%[~, NumRegionsOriginalOtsuMask] = bwlabel(OriginalOtsuMask, 8);
% do a second LK on the probability map
%%


%%
%TemplateOverlapOtsu = sum(OtsuMask(:) .* PIMGCropped(:)) ./ sum(PIMGCropped(:));
%[~, OriginalOtsuSeg] = histc(ResampledAVWCroppedForOtsu(:), [0, T, 256]);
%OriginalOtsuSeg = reshape(OriginalOtsuSeg, size(ResampledAVWCroppedForOtsu));
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

OtsuMaskCC = bwconncomp(OtsuMask);
L = labelmatrix(OtsuMaskCC);
R = regionprops(OtsuMaskCC, 'Area', 'PixelIdxList');

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

OtsuMask = ismember(L, find(TemplateOverlap > 0.25));
clear Junk L;

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

FirstStartMask = OtsuMask & (PIMGCropped > 0.5);

MeanIntensityInStartMask = mean(ResampledAVWCropped(FirstStartMask));
VarIntensityInStartMask = var(ResampledAVWCropped(FirstStartMask));

GaussianProb = normpdf(ResampledAVWCropped, MeanIntensityInStartMask, sqrt(VarIntensityInStartMask));
GaussianProb = GaussianProb ./ max(GaussianProb(:));
GaussianProb(ResampledAVWCropped > MeanIntensityInStartMask) = 1;

SIGMA = 0.25;

[GaussianFilter, GaussianFilterDeriv] = gaussian_filter_max_1d(SIGMA);

aSmooth=imfilter(GaussianProb, GaussianFilter, 'conv', 'replicate', 'same');   % run the filter across rows
aSmooth=imfilter(aSmooth, GaussianFilter','conv','replicate', 'same'); % and then across columns

%apply directional derivatives
ax = imfilter(aSmooth, GaussianFilterDeriv(:)', 'conv','replicate');
ay = imfilter(aSmooth, GaussianFilterDeriv(:), 'conv','replicate');

GaussianGradMAG = sqrt(ax .* ax + ay .* ay);
GaussianGradMAG = GaussianGradMAG ./ max(GaussianGradMAG(:));
%ExclusionWeight = 
ExclusionA = min(FIMGCropped * 50, 1) .* GaussianProb;
ExclusionB = GaussianGradMAG;
ExclusionC = exp(-2 .* GaussianProb .* (PIMGCropped + 2));
ExclusionA = ExclusionA ./ max(ExclusionA(:));
ExclusionB = ExclusionB ./ max(ExclusionB(:));
ExclusionC = ExclusionC ./ max(ExclusionC(:));

%ExclusionWeight = min(FIMGCropped * 50, 1) .* GaussianProb + GaussianGradMAG + exp(-2*GaussianProb);
%ExclusionWeight = ExclusionWeight ./ max(ExclusionWeight(:));
ExclusionWeight = ExclusionA + 4 * ExclusionB + 2 * ExclusionC;
%StartMask = StartMask & ExclusionWeight < 0.25;
%keyboard;
% subplot 324;
% imshow(P, []);
% subplot 326;
% use hysteresis on the exclusion mask to get an initial case
% strong is the low exclusion weight that overlaps with the template
% weak is just the low exclusion weight
% find a new start mask with hysteresis

% clf;
% %[StrongI, StrongJ] = 
% %subplot 121;
% imshow(ExclusionWeight, []);
% hold on;
% [C, CC] = contour(FirstStartMask, [0.5 0.5]);
% set(CC, 'Color', 'r');

[StrongI, StrongJ] = find(FirstStartMask);
StartMask = bwselect(ExclusionWeight < 0.8, StrongJ, StrongI, 8);

%L = regionprops(StartMask, 'Area');
%StartMaskAreas = [L.Area];

StartMask = bwareaopen(StartMask, 200);
%subplot 122;
% imshow(StartMask, []);
OriginalPSI = double(StartMask);
%OriginalS = OtsuMask & (PIMGCropped > 0.2) & (GaussianProb > 0.3 | ResampledAVWCropped > MeanIntensityInStartMask);
%StartMask;

%T = S + FIMGCropped;
%imshow(T, []);
% subplot 321;
% hold on;
% [L, CC] = contour(S, [0.5, 0.5]);
%keyboard;

DMask = (GaussianProb > 0.2 | ResampledAVWCropped > MeanIntensityInStartMask) & ExclusionWeight < 0.5;
DMask = imdilate(DMask, strel('disk', 3));

% [L, CCT] = contour(DMask, [0.5, 0.5]);
% set(CCT, 'Color', 'r');


StartMask = bwfill(StartMask, 'holes', 8);
StartMaskCC = bwconncomp(StartMask, 8);
NumRegionsStartMask = StartMaskCC.NumObjects;
%[L, NumRegionsStartMask] = bwlabel(StartMask, 8);

if(NumRegionsStartMask > 1)
	OldStartMask = StartMask;
	
	[B] = bwboundaries(StartMask, 8);
	Distances = zeros(length(B));
	DistancesMinI = zeros(length(B));
	DistancesMinJ = zeros(length(B));
	
	for BoundaryI = 1:length(B) - 1
		for BoundaryJ = BoundaryI + 1:length(B)
			XC = bsxfun(@minus, B{BoundaryI}(:, 1), B{BoundaryJ}(:, 1)');
			YC = bsxfun(@minus, B{BoundaryI}(:, 2), B{BoundaryJ}(:, 2)');
			SZ = size(XC);
			XC = XC(:);
			YC = YC(:);
			[Distances(BoundaryI, BoundaryJ), I] = min(sqrt(XC .* XC + YC .* YC));
			[DistancesMinI(BoundaryI, BoundaryJ), DistancesMinJ(BoundaryI, BoundaryJ)] = ind2sub(SZ, I);
		end
	end
	Distances = Distances + Distances';
	DistancesMinI = DistancesMinI + DistancesMinI';
	DistancesMinJ = DistancesMinJ + DistancesMinJ';
	if(length(B) > 2)
		edges = prim(Distances);
	else
		edges = [1, 2];
	end

	for z = 1:size(edges, 1)
		StartX = B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 2);
		EndX = B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 2);
		StartY = B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 1);
		EndY = B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 1);
		
		ArcLength = sqrt((StartX - EndX) * (StartX - EndX) + (StartY - EndY) * (StartY - EndY));
		n = ceil(ArcLength * sqrt(2)); % make sure we cover all pixels
		
		IX = linspace(StartX, EndX, n);
		IY = linspace(StartY, EndY, n);
		IX = round(IX);
		IY = round(IY);
		I = sub2ind(size(OldStartMask), IY, IX);
		I = unique(I);
		T = false(size(StartMask));
		T(I) = 1;
		%clear I IX IY n ArcLength EndY StartY EndX StartX;
% 		line([...
% 		B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 2), ...
% 		B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 2)], ... 
% 		[B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 1), ...
% 		B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 1)]);
% 	
		Angle = atan2(EndY - StartY, EndX - StartX);
		Angle = mod(Angle + 2 * pi, 2 * pi);

		R = [5, 5] + 7 * abs([cos(Angle), sin(Angle)]);

		%[GX(I, J), GY(I, J)]

		AngleWeighting = -0.9 * cos(2 * (Angle + 45 * pi / 180));

		SQ = sqrt(R(1) * R(2)) * AngleWeighting;
		SIGMA = [R(1), SQ; SQ, R(2)];
		[F, HalfMaximum] = gaussian_fwhm2d(SIGMA);
		
		T = imdilate(T, double(F > HalfMaximum));
		StartMask = StartMask | (T & GaussianProb > 0.3);
		DMask = DMask | (imdilate(T, strel('disk', 2)) & GaussianProb > 0.2);
% 		subplot 131;
% 		imshow(OldStartMask, []);
% 		subplot 132;
% 		imshow(StartMask, []);
% 		subplot 133;
% 		imshow(F > HalfMaximum, []);
% 		keyboard;
		
		clear I IX IY n ArcLength EndY StartY EndX StartX F T HalfMaximum SQ Angle AngleWeighting R;
	end
% 	subplot 121;
% 	imshow(OldStartMask, []);
% 	subplot 122;
% 	imshow(StartMask, []);
% 	keyboard;
	UsedConvexHull = true;
% 	return;
	%[D, C] = radial_gradient_dilate_fwhm(StartMask, double(AIMGCropped), false);
	%[DE, C] = radial_gradient_dilate_fwhm(~D, double(AIMGCropped), true);
	%DE = ~DE;
	%%
% 	clf;
% 	XX = 65;
% 	YY = 42;
% % 	subplot 121;
% % 	imshow(StartMask, []);
% % 	hold on;
% % 	plot(XX, YY, '*');
% 	%[XX, YY] = getpts;
% 	%XX = round(XX);
% 	%YY = round(YY);
% 	%subplot 122;
% 	E = exp(1i * AIMGCropped);
% 	subplot 121;
% 	hold off;
% 	imshow(StartMask, []);
% 	hold on;
% 	[X, Y] = meshgrid(1:size(E, 2), 1:size(E, 1));
% 	quiver(X, Y, real(E), imag(E), 0);
% 	B = bwboundaries(StartMask);
% 	
% 	SForward = cell(1, length(B));
% 	SBackward = cell(1, length(B));
% 	S = cell(1, length(B));
% 	for z = 1:length(B)
% 		SForward{z} = stream2(X, Y, real(E), imag(E), B{z}(:, 2), B{z}(:, 1), [1, 50]);
% 		SBackward{z} = stream2(X, Y, -real(E), -imag(E), B{z}(:, 2), B{z}(:, 1), [1, 50]);
% 		S{z} = cell(length(SForward{z}), 1);
% 		for k = 1:length(SForward{z})
% 			S{z}{k} = [SForward{z}{k}(end:-1:2, :); SBackward{z}{k}];
% 		end
% 	end
% 	axis equal ij;
% 	%zoom(2.5);
% 	%hold on;
% 	subplot 122;
% 	imshow(StartMask, []);
% 	Colors = 'rb';
% 	for z = 1:length(B)
% 		L = streamline(S{z});
% 		set(L, 'Color', Colors(z));
% 	end
% 	
% 	%plot(XX, YY, 'r*');
% 	%imshow(AIMGCropped, []);
% % 	subplot 223;
% % 	imshow(D, []);
% % 	subplot 224;
% % 	imshow(DE, []);
% 	%%
% 	keyboard;
	%%
%	disp('Using convex hull');
% 	[Y, X] = find(StartMask);
% 	[ConvHullX] = convhull(X, Y);
% 	BW = roipoly(StartMask, X(ConvHullX), Y(ConvHullX));
% 	
% 	OriginalPSI = zeros(size(BW));
% 	OriginalPSI(BW) = 1;
% 	OriginalPSI(~BW) = -1;
% 	
% 	OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);
% 	
% 	D = double(StartMask);
% 	D(StartMask) = 1;
% 	D(~StartMask) = -1;
% 	
% 	D = computeDistanceFunction2d(D, [1, 1], [], 2);
% 	
% 	[NewPSI, NumIter] = fls_evolution_illusory_contours(-D, OriginalPSI, ones(size(BW)), zeros(size(BW)), [1, 1], 1000, 200, 5, 50, ResampledAVWCropped, true);
% 	keyboard;
% 	BW = imdilate(BW, strel('disk', 4));
% 	OldStartMask = StartMask;
% 	StartMask = StartMask | (BW & PIMGCropped > 0.15);
% 	StartMask = bwfill(StartMask, 'holes', 8);
% 	AddedByConvexHull = StartMask & ~OldStartMask;
	
	
	
% 	disp('Using dilation method');
% 	
% 	DilationR = [3, 7];
% 	%DilationR = 3;
% 	while(1)
% 		T = padarray(OldStartMask, DilationR, 0, 'both');
% 		
% 		DilationStartMask = radial_close(T, DilationR);
% 		%T = radial_dilate(T, DilationR - 3);
% 		%SE = strel('disk', DilationR);
% 		
% 		DilationStartMask = DilationStartMask(DilationR(1) + 1:DilationR(1) + size(OldStartMask, 1), DilationR(2) + 1:DilationR(2) + size(OldStartMask, 2));
% 		%DilationStartMask = imclose(OldStartMask, SE);
% 		DilationStartMask = DilationStartMask | OldStartMask;
% 		[~, NumLabels] = bwlabel(DilationStartMask, 8);
% 		imshow(DilationStartMask, []);
% 		title(num2str(DilationR));
% 		keyboard;
% 		if(NumLabels == 1)
% 			disp(['Dilation size: ' num2str(DilationR)]);
% 			break;
% 		end
% 
% 		DilationR = DilationR + [2, 4];
% 		
% 	end
else
	UsedConvexHull = false;
	%return;
end

% clf;
% if(UsedConvexHull)
% 	subplot 221;
% 	imshow(ResampledAVWCropped, []);
% 	hold on;
% 	[C, CC] = contour(StartMask, [0.5 0.5]);
% 	set(CC, 'Color', 'r');
% 	[C, CC] = contour(OldStartMask, [0.5 0.5]);
% 	set(CC, 'Color', 'g');
% 	%title(num2str(StartMaskAreas));
% 	subplot 222;
% 	imshow(PIMGCropped, []);
% 	hold on;
% 	[C, CC] = contour(StartMask, [0.5 0.5]);
% 	set(CC, 'Color', 'r');
% 	[C, CC] = contour(OldStartMask, [0.5 0.5]);
% 	set(CC, 'Color', 'g');
% 	[C, CC] = contour(PIMGCropped, [0.5 0.5]);
% 	set(CC, 'Color', 'b');
% 	
% % 	subplot 223;
% % 	imshow(AddedByConvexHull, []);
% 	subplot 224;
% 	imshow(ResampledAVWCropped, []);
% 	hold on;
% % 	[C, CC] = contour(DilationStartMask, [0.5, 0.5]);
% % 	set(CC, 'Color', 'r');
% 	[C, CC] = contour(OldStartMask, [0.5, 0.5]);
% 	set(CC, 'Color', 'b');
% else
% 	imshow(ResampledAVWCropped, []);
% 	hold on;
% 	[C, CC] = contour(StartMask, [0.5 0.5]);
% 	set(CC, 'Color', 'b');
% end
% 
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

% return;
% 	
% if(NumRegionsStartMask > 1)
% 	SR = 2;
% 	SC = 3;
% 	subplot(SR, SC, 1);
% 	hold off;
% 	imshow(L, []);
% 	[Y, X] = find(StartMask);
% 	[ConvHullX] = convhull(X, Y);
% 	hold on;
% 	plot(X(ConvHullX), Y(ConvHullX));
% 	title('L');
% 	subplot(SR, SC, 2);
% 	imshow(DMask, []);
% 	title('DMask');
% 	subplot(SR, SC, 3);
% 	%imshow(ExclusionWeight, []);
% 	imshow(ResampledAVWCropped, []);
% 	%R = regionprops(StartMask, 'ConvexImage');
% 	%keyboard;
% 	%imshow(R.ConvexImage, []);
% %	C
% 	title('OtsuMask');
% 	subplot(SR, SC, 4);
% 	%%
% 	BW = roipoly(StartMask, X(ConvHullX), Y(ConvHullX));
% 	BW = imdilate(BW, strel('disk', 4));
% 	StartMask = StartMask | (BW & PIMGCropped > 0.5);
% 	%StartMask = StartMask & (ExclusionB < 0.5);
% 	imshow(StartMask, []);
% 	title('New StartMask');
% 	subplot(SR, SC, 3);
% 	%imshow(ExclusionWeight, []);
% 	imshow(ResampledAVWCropped, []);
% 	hold on;
% 	[~, CC] = contour(StartMask, [0.5 0.5]);
% 	
% 	set(CC, 'Color', 'b');
% 	%%
% 	subplot(SR, SC, 5);
% 	
% 	imshow(ExclusionWeight, []);
% 	hold on;
% 	[~, CC] = contour(DMask, [0.5 0.5]);
% 	
% 	set(CC, 'Color', 'b');
% 	subplot(SR, SC, 6);
% 	imshow(ExclusionC, []);
% 	
% 	return;
% end
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
if(all(OriginalPSI(:) < 0) || all(OriginalPSI(:) > 0))
%if(all(D(:) > 0) || all(D(:) < 0) || any(isnan(GaussianProb(:))))
	ContourLevel = 0.15;
	clf;
	subplot 221;
	imshow(ResampledAVW, []);
	hold on;
	SZ = size(MatchedFilter);

	T = line([realxpeak, realxpeak + SZ(2)], [realypeak, realypeak]); set(T, 'Color', 'b');
	T = line([realxpeak, realxpeak + SZ(2)], [realypeak + SZ(1), realypeak + SZ(1)]); set(T, 'Color', 'b');
	T = line([realxpeak, realxpeak], [realypeak, realypeak + SZ(1)]); set(T, 'Color', 'b');
	T = line([realxpeak + SZ(2), realxpeak + SZ(2)], [realypeak, realypeak + SZ(1)]); set(T, 'Color', 'b');

	%plot(TX(:, 1), TY(:, 1), 'r*');
	%plot(TX(:, end), TY(:, end), 'r*');
	%plot(TX(1, :), TY(1, :), 'r*');
	%plot(TX(end, :), TY(end, :), 'r*');
	
	%plot(TX, TY, '*');
	title('Initial XCORR alignment');
	subplot 222;
	imshow(ResampledAVW, []);
	hold on;
	SZ = size(MatchedFilter);
	
	T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
	T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
	T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
	T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
	
	[~, CC] = contour(PIMG, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
	%[~, CC] = contour(PIMGW, [ContourLevel, ContourLevel]); set(CC, 'Color', 'y');
	
	title('Initial LK alignment, red uniform, yellow weighted');
	keyboard;
% 	subplot 223;
% 	imshow(ResampledAVW, []);
% 	hold on;
% 	SZ = size(MatchedFilter);
% 	
% 	T = line([TXSeg(1, 1) TXSeg(end, 1)], [TYSeg(1, 1), TYSeg(end, 1)]); set(T, 'Color', 'r');
% 	T = line([TXSeg(1, end) TXSeg(end, end)], [TYSeg(1, end), TYSeg(end, end)]); set(T, 'Color', 'r');
% 	T = line([TXSeg(1, 1) TXSeg(1, end)], [TYSeg(1, 1), TYSeg(1, end)]); set(T, 'Color', 'r');
% 	T = line([TXSeg(end, 1) TXSeg(end, end)], [TYSeg(end, 1), TYSeg(end, end)]); set(T, 'Color', 'r');
% 	
% 	T = line([TXSegW(1, 1) TXSegW(end, 1)], [TYSegW(1, 1), TYSegW(end, 1)]); set(T, 'Color', 'y');
% 	T = line([TXSegW(1, end) TXSegW(end, end)], [TYSegW(1, end), TYSegW(end, end)]); set(T, 'Color', 'y');
% 	T = line([TXSegW(1, 1) TXSegW(1, end)], [TYSegW(1, 1), TYSegW(1, end)]); set(T, 'Color', 'y');
% 	T = line([TXSegW(end, 1) TXSegW(end, end)], [TYSegW(end, 1), TYSegW(end, end)]); set(T, 'Color', 'y');
% 
% 	[~, CC] = contour(PIMGSeg, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
% 	[~, CC] = contour(PIMGSegW, [ContourLevel, ContourLevel]); set(CC, 'Color', 'y');
% 	title('Initial LK alignment to seg, red uniform, yellow weighted');
% 	
% % 	subplot 223;
% % 	%RWeightingIMG = padarray(RWeighting, size(MatchedFilter), 0, 'pre');
% % 	imshow(cc2R, []);
% % 	title('Correlation coefficient output');
% 	subplot 224;
% 	%RWeightingIMG = padarray(RWeighting, size(MatchedFilter), 0, 'pre');
% 	imshow(PIMG, []);
% 	hold on;
% 	T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
% 	T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
% 	T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
% 	T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
% 	
% 	[~, CC] = contour(PIMG, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
% 	
	%imshow(WMSegSmoothed, []);
	%title('WMSeg Smoothed');
	title('P');
	FigPos = fullscreen_fig_pos;
	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
	exportfig(gcf, [OutputPNG '-cc2'], 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

	disp('Bad image, returning');
	return;
	
end
%if(all(D(:) > 0) || all(D(:) < 0) || any(isnan(GaussianProb(:))))
% ContourLevel = 0.15;
% 	clf;
% 	subplot 221;
% 	imshow(ResampledAVW, []);
% 	hold on;
% 	SZ = size(MatchedFilter);
% 
% 	T = line([xpeak, xpeak + SZ(2)], [ypeak, ypeak]); set(T, 'Color', 'b');
% 	T = line([xpeak, xpeak + SZ(2)], [ypeak + SZ(1), ypeak + SZ(1)]); set(T, 'Color', 'b');
% 	T = line([xpeak, xpeak], [ypeak, ypeak + SZ(1)]); set(T, 'Color', 'b');
% 	T = line([xpeak + SZ(2), xpeak + SZ(2)], [ypeak, ypeak + SZ(1)]); set(T, 'Color', 'b');
% 
% 	T = line([xpeak2, xpeak2 + SZ(2)], [ypeak2, ypeak2]); set(T, 'Color', 'g');
% 	T = line([xpeak2, xpeak2 + SZ(2)], [ypeak2 + SZ(1), ypeak2 + SZ(1)]); set(T, 'Color', 'g');
% 	T = line([xpeak2, xpeak2], [ypeak2, ypeak2 + SZ(1)]); set(T, 'Color', 'g');
% 	T = line([xpeak2 + SZ(2), xpeak2 + SZ(2)], [ypeak2, ypeak2 + SZ(1)]); set(T, 'Color', 'g');
% 
% 	T = line([xpeak2R, xpeak2R + SZ(2)], [ypeak2R, ypeak2R]); set(T, 'Color', 'm');
% 	T = line([xpeak2R, xpeak2R + SZ(2)], [ypeak2R + SZ(1), ypeak2R + SZ(1)]); set(T, 'Color', 'm');
% 	T = line([xpeak2R, xpeak2R], [ypeak2R, ypeak2R + SZ(1)]); set(T, 'Color', 'm');
% 	T = line([xpeak2R + SZ(2), xpeak2R + SZ(2)], [ypeak2R, ypeak2R + SZ(1)]); set(T, 'Color', 'm');
% 
% 	%plot(TX(:, 1), TY(:, 1), 'r*');
% 	%plot(TX(:, end), TY(:, end), 'r*');
% 	%plot(TX(1, :), TY(1, :), 'r*');
% 	%plot(TX(end, :), TY(end, :), 'r*');
% 	
% 	%plot(TX, TY, '*');
% 	title('Initial XCORR alignment, blue entire image, green seg, magenta seg with weight');
% 	subplot 222;
% 	imshow(ResampledAVW, []);
% 	hold on;
% 	SZ = size(MatchedFilter);
% 	
% 	T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
% 	T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
% 	T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
% 	T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
% 	
% 	T = line([TXW(1, 1) TXW(end, 1)], [TYW(1, 1), TYW(end, 1)]); set(T, 'Color', 'y');
% 	T = line([TXW(1, end) TXW(end, end)], [TYW(1, end), TYW(end, end)]); set(T, 'Color', 'y');
% 	T = line([TXW(1, 1) TXW(1, end)], [TYW(1, 1), TYW(1, end)]); set(T, 'Color', 'y');
% 	T = line([TXW(end, 1) TXW(end, end)], [TYW(end, 1), TYW(end, end)]); set(T, 'Color', 'y');
% 	[~, CC] = contour(PIMG, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
% 	[~, CC] = contour(PIMGW, [ContourLevel, ContourLevel]); set(CC, 'Color', 'y');
% 	
% 	title('Initial LK alignment, red uniform, yellow weighted');
% 	
% 	subplot 223;
% 	imshow(ResampledAVW, []);
% 	hold on;
% 	SZ = size(MatchedFilter);
% 	
% 	T = line([TXSeg(1, 1) TXSeg(end, 1)], [TYSeg(1, 1), TYSeg(end, 1)]); set(T, 'Color', 'r');
% 	T = line([TXSeg(1, end) TXSeg(end, end)], [TYSeg(1, end), TYSeg(end, end)]); set(T, 'Color', 'r');
% 	T = line([TXSeg(1, 1) TXSeg(1, end)], [TYSeg(1, 1), TYSeg(1, end)]); set(T, 'Color', 'r');
% 	T = line([TXSeg(end, 1) TXSeg(end, end)], [TYSeg(end, 1), TYSeg(end, end)]); set(T, 'Color', 'r');
% 	
% 	T = line([TXSegW(1, 1) TXSegW(end, 1)], [TYSegW(1, 1), TYSegW(end, 1)]); set(T, 'Color', 'y');
% 	T = line([TXSegW(1, end) TXSegW(end, end)], [TYSegW(1, end), TYSegW(end, end)]); set(T, 'Color', 'y');
% 	T = line([TXSegW(1, 1) TXSegW(1, end)], [TYSegW(1, 1), TYSegW(1, end)]); set(T, 'Color', 'y');
% 	T = line([TXSegW(end, 1) TXSegW(end, end)], [TYSegW(end, 1), TYSegW(end, end)]); set(T, 'Color', 'y');
% 
% 	[~, CC] = contour(PIMGSeg, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
% 	[~, CC] = contour(PIMGSegW, [ContourLevel, ContourLevel]); set(CC, 'Color', 'y');
% 	title('Initial LK alignment to seg, red uniform, yellow weighted');
% 	
% % 	subplot 223;
% % 	%RWeightingIMG = padarray(RWeighting, size(MatchedFilter), 0, 'pre');
% % 	imshow(cc2R, []);
% % 	title('Correlation coefficient output');
% 	subplot 224;
% 	%RWeightingIMG = padarray(RWeighting, size(MatchedFilter), 0, 'pre');
% 	imshow(PIMG, []);
% 	hold on;
% 	T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
% 	T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
% 	T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
% 	T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
% 	
% 	[~, CC] = contour(PIMG, [ContourLevel, ContourLevel]); set(CC, 'Color', 'r');
% 	
% 	%imshow(WMSegSmoothed, []);
% 	%title('WMSeg Smoothed');
% 	title('P');
% 	FigPos = fullscreen_fig_pos;
% 	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% 	exportfig(gcf, [OutputPNG '-cc2'], 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 	disp('Bad image, returning');
% 	return;
% %else
% %	return;
% %end
% return;

% if(UsedConvexHull)
% % 	clf;
% % 	subplot 221;
% % 	imshow(ResampledAVWCropped, []);
% % 	hold on;
% % 	[~, CC] = contour(DMask, [0.5, 0.5]);
% % 	set(CC, 'Color', 'r');
% % 	[~, CC] = contour(StartMask, [0.5, 0.5]);
% % 	set(CC, 'Color', 'g');
% % 	subplot 222;
%  	
% % 	subplot 223;
% 	%DMaskL = bwlabel(DMask);
% 	%imshow(DMask, []);
% 	
% 	
% 	
% % 	subplot 224;
% % 	imshow(ResampledAVWCropped, []);
% % 	hold on;
% % 	[~, CC] = contour(DMask, [0.5, 0.5]);
% % 	set(CC, 'Color', 'r');
% % 	[~, CC] = contour(StartMask, [0.5, 0.5]);
% % 	set(CC, 'Color', 'g');
% else
% 	delete(OutputPNG);
% end

DMaskCC = bwconncomp(DMask);
DMaskL = labelmatrix(DMaskCC);
% 	imshow(label2rgb(DMaskL, 'lines', 'k'), []);
% 	
% find the region in DMask that mostly overlaps with the startmask, filter out other regions
DMaskR = regionprops(DMaskL, 'PixelIdxList');

Overlaps = zeros(1, length(DMaskR));
I = find(StartMask);
for z = 1:length(DMaskR)
	TF = ismember(DMaskR(z).PixelIdxList, I);
	Overlaps(z) = sum(TF) ./ length(I);
	clear TF;
end
clear I;
if(all(Overlaps < 0.975))
	clf;
	imshow(label2rgb(DMaskL, 'lines', 'k'), []);
	keyboard;
else
	[~, L] = max(Overlaps);
	DMask = (DMaskL == L);
end
clear Overlaps DMaskL L;

SE = strel('disk', 4);
DMask = imclose(DMask, SE);
StartMask = imclose(StartMask, SE);

OriginalPSI = zeros(size(StartMask));
OriginalPSI(StartMask) = 1;
OriginalPSI(~StartMask) = -1;
%OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);
OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);

% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% return;

D = double(DMask);
D(DMask) = 1;
D(~DMask) = -1;
D = computeDistanceFunction2d(D, [1, 1], [], 2);
%clear S;

%keyboard;
%FIMGCropped(isinf(FIMGCropped)) = 0;
W = 5;
%[NewPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(-D, OriginalPSI, W - (PIMGCropped * (W - 1)), double(ExclusionWeight), [1, 1], 1000, 200, 10, 5, ResampledAVWCropped);
%[NewPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(-D, OriginalPSI, ones(size(D)), double(ExclusionWeight), [1, 1], 1000, 200, 10, 5, ResampledAVWCropped);
[NewPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(-D, OriginalPSI, ones(size(D)), zeros(size(D)), [1, 1], 1000, 200, 10, 5, ResampledAVWCropped);
%keyboard;
PlotGraphics = true;

if(PlotGraphics)
% 	clf;
% 	SR = 3;
% 	SC = 3;
% 	%keyboard;
% 	subplot(SR, SC, 1);
% 	imshow(cc2R, []);
% 	hold on;
% 	plot(SortedCCJD(FirstOffsetValid2R), SortedCCID(FirstOffsetValid2R), 'r*');
% 	plot(cc2RRegMaxIMaxJ, cc2RRegMaxIMaxI, 'b*');
% 
% 	title('Normxcorr result, red based on maximum cc, blue based on minimum distance from expected');
% 
% 	subplot(SR, SC, 2);
% 	imshow(cc2R, []);
% 	%[I, J] = find
% 	hold on;
% 	plot(cc2RRegMaxJ, cc2RRegMaxI, '*');
% 	title('Regional maxima of cc2R');
% 
% 	subplot(SR, SC, 3);
% 	imshow(ResampledAVW, []);
% 	hold on;
% 	SZ = size(MatchedFilter);
% 
% 	T = line([realxpeak, realxpeak + SZ(2)], [realypeak, realypeak]); set(T, 'Color', 'b');
% 	T = line([realxpeak, realxpeak + SZ(2)], [realypeak + SZ(1), realypeak + SZ(1)]); set(T, 'Color', 'b');
% 	T = line([realxpeak, realxpeak], [realypeak, realypeak + SZ(1)]); set(T, 'Color', 'b');
% 	T = line([realxpeak + SZ(2), realxpeak + SZ(2)], [realypeak, realypeak + SZ(1)]); set(T, 'Color', 'b');
% 	% T = line([realxpeak2, realxpeak2 + SZ(2)], [realypeak2, realypeak2]); set(T, 'Color', 'g');
% 	% T = line([realxpeak2, realxpeak2 + SZ(2)], [realypeak2 + SZ(1), realypeak2 + SZ(1)]); set(T, 'Color', 'g');
% 	% T = line([realxpeak2, realxpeak2], [realypeak2, realypeak2 + SZ(1)]); set(T, 'Color', 'g');
% 	% T = line([realxpeak2 + SZ(2), realxpeak2 + SZ(2)], [realypeak2, realypeak2 + SZ(1)]); set(T, 'Color', 'g');
% 
% 	%plot(TX(:, 1), TYSeg(:, 1), 'r*');
% 	%plot(TX(:, end), TY(:, end), 'r*');
% 	%plot(TX(1, :), TY(1, :), 'r*');
% 	%plot(TX(end, :), TY(end, :), 'r*');
% 	T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
% 	T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
% 	T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
% 	T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
% 	%plot(TX, TY, '*');
% 	title('Initial XCORR and LK alignment');
% 
% 	subplot(SR, SC, 4);
% 	imshow(ResampledAVWCropped, []);
% 	hold on;
% 	[~, CC] = contour(NewPSI, [0, 0]);
% 	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
% 	title('Original Image');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
% 
% 	subplot(SR, SC, 5);
% 	imshow(GaussianProb, []);
% 	hold on;
% 	[~, CC] = contour(NewPSI, [0, 0]);
% 	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
% 	title('Gaussian probability of pixels');
% 
% 	subplot(SR, SC, 6);
% 	imshow(PIMGCropped, []);
% 	hold on;
% 	[~, CC] = contour(NewPSI, [0, 0]);
% 	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
% 	title('Probability of CC from template');
% 
% 	subplot(SR, SC, 7);
% 	imshow(PIMGCropped, []);
% 	hold on;
% 	[~, CC] = contour(OriginalOtsuMask, [0.5, 0.5]);
% 	set(CC, 'Color', 'r');
% 	%title(['Probability of CC from template with original otsu borders, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']);
% 	title(['Probability of CC from template with original otsu borders']);
% 
% % 	subplot(SR, SC, 8);
% % 	imshow(OriginalOtsuSeg, []);
% % 	hold on;
% % 	[~, CC] = contour(NewPSI, [0, 0]);
% % 	set(CC, 'Color', 'r');
% % 	[~, CC] = contour(OriginalPSI, [0, 0]);
% % 	set(CC, 'Color', 'b');
% % 	title(['OtsuMask ' num2str(OriginalOtsuSegCounts(1:3)', '%d ')]);
% 
% 	subplot(SR, SC, 9);
% 	imshow(ExclusionWeight, []);
% 	hold on;
% 	[~, CC] = contour(NewPSI, [0, 0]);
% 	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
% 	title('Exclusion mask');

	clf;
	
	imshow(ResampledAVWCropped, []);
	hold on;
	[~, CC] = contour(NewPSI, [0, 0]);
	set(CC, 'Color', 'r');
	[~, CC] = contour(OriginalPSI, [0, 0]);
	set(CC, 'Color', 'b');
	title('Original Image');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']

	FigPos = fullscreen_fig_pos;
	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
	exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
	
	IMG = imread(OutputPNG);
	IMG = imautocropwhite(IMG);
	imwrite(IMG, OutputPNG);
end

OutputS.InitialSeg = OriginalPSI > 0;
OutputS.FinalSeg = NewPSI > 0;
OutputS.GroundSeg = ResampledGroundCropped > 0;
OutputS.InitialDice = dices_coefficient(OutputS.InitialSeg, OutputS.GroundSeg);
OutputS.FinalDice = dices_coefficient(OutputS.FinalSeg, OutputS.GroundSeg);
[~, OutputS.NumRegions] = bwlabel(OutputS.FinalSeg);

save(OutputMAT, '-struct', 'OutputS');

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


function [D] = dices_coefficient(A, B)

D = 2 .* sum(A(:) & B(:)) ./ (sum(A(:)) + sum(B(:)));
