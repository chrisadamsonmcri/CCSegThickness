function [ReturnCode] = cc_seg_one_subject_pipe_seg_python(MATFilePrefix, GroundTruthFile, DoLK)

% cc_seg_one_subject_pipe_seg(MATFilePrefix, GroundTruthFile)
%
% DESCRIPTION
%	Performs corpus callosum segmentation the using illusory contours level set method and various mathematical morphology operations.
%	Performs the segmentation part only.
%
% PARAMETERS
% MATFilePrefix (string) [1]: the
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%

[pathstr, name, ext] = fileparts(MATFilePrefix);

OutputDir = pathstr;
OutputPrefix = name;
if(strcmp(ext, '.mat'))
	FileName = name;
else
	FileName = [name ext];
end

if(exist(fullfile(pathstr, [FileName '_midsag.mat']), 'file') == 2)
	MidsagMatFileName = fullfile(pathstr, [FileName '_midsag.mat']);
	MidsagMatFileExt = 'mat';
elseif(exist(fullfile(pathstr, [FileName '_midsag.hdf5']), 'file') == 2)
	MidsagMatFileName = fullfile(pathstr, [FileName '_midsag.hdf5']);
	MidsagMatFileExt = 'hdf5';
else
	error('Midsag file not found');
end

switch(MidsagMatFileExt)
	case 'mat'
		load(MidsagMatFileName, 'AVWMidSag', 'NIIPixdims');
	case 'hdf5'
		%keyboard;
		
		% newer matlabs use h5read, try this first
		try
			AVWMidSag = h5read(MidsagMatFileName, '/midSagAVW');
			NIIPixdims = h5read(MidsagMatFileName, '/NIIPixdims');
		catch e
			% otherwise, use the older hdf5read
			%I = hdf5info(MidsagMatFileName);
			AVWMidSag = hdf5read(MidsagMatFileName, '/midSagAVW');
			NIIPixdims = hdf5read(MidsagMatFileName, '/NIIPixdims');
			%keyboard;
		end
		AVWMidSag = flipdim(AVWMidSag, 2);
		AVWMidSag = rot90(AVWMidSag, 1);
		%whos AVWMidSag
		%keyboard;
end

SegMatFileName = fullfile(pathstr, [FileName '_seg.mat']);

GroundTruthGiven = false;
if(~isempty(GroundTruthFile) && ischar(GroundTruthFile))
	if(exist(GroundTruthFile, 'file') == 2)
		[~, GroundTruthAVW] = load_nii(GroundTruthFile);
		GroundTruthGiven = true;
		disp('Using ground truth for comparison');
	else
		error('Ground truth file was given but does not exist');
	end
end

% load atlases

if(exist('cc_seg_atlases.mat', 'file') == 2)
	load('cc_seg_atlases.mat');
else
	disp('Loading atlases from NIFTIs');
	AtlasDir = '/data/addo/cc_seg';
	[~, FornixProbAVW] = load_nii(fullfile(AtlasDir, 'all_fornix_prob.nii.gz'));
	[~, ProbAVW] = load_nii(fullfile(AtlasDir, 'all_cc_prob.nii.gz'));
	[TemplateNII, TemplateAVW] = load_nii(fullfile(AtlasDir, 'all_msp_mean.nii.gz'));
	
	save('cc_seg_atlases.mat', 'FornixProbAVW', 'ProbAVW', 'TemplateNII', 'TemplateAVW');
end
TemplateAVW = double(TemplateAVW);
ProbAVW = double(ProbAVW);
FornixProbAVW = double(FornixProbAVW);

TemplatePixdims = TemplateNII.hdr.dime.pixdim(3:4);

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

MatchedFilter = imcrop(TemplateAVW, R);
MatchedFilterProb = imcrop(ProbAVW, R);
MatchedFilterFornixProb = imcrop(FornixProbAVW, R);
%keyboard;
LeftLimit = max(min(Templatexx), min(AVWMidSagxx));
UpperLimit = max(min(Templateyy), min(AVWMidSagyy));
RightLimit = min(max(Templatexx), max(AVWMidSagxx));
LowerLimit = min(max(Templateyy), max(AVWMidSagyy));

[ResampleX, ResampleY] = meshgrid(LeftLimit:TemplatePixdims(1):RightLimit, UpperLimit:TemplatePixdims(2):LowerLimit);

ResampledAVW = interp2(AVWMidSagX, AVWMidSagY, AVWMidSag, ResampleX, ResampleY);
ResampledAVW = double(ResampledAVW);
if(GroundTruthGiven)
	ResampledGroundAVW = interp2(AVWMidSagX, AVWMidSagY, GroundTruthAVW, ResampleX, ResampleY, 'nearest');
else
	ResampledGroundAVW = [];
end

[THRESH, SEG] = robust_otsu2(ResampledAVW(ResampledAVW > 0), [0.05, 0.98]);
if(isempty(THRESH))
	disp('acpcdetect failed');
	ReturnCode = 1;
	return;
end
FullSEG = zeros(size(ResampledAVW));
FullSEG(ResampledAVW > 0) = SEG;
FullSEG(ResampledAVW > THRESH(end)) = 3;
%FullOtsuSeg = FullSEG;
WMSeg = (FullSEG == 3);
%keyboard;
[I, J] = find(WMSeg);
clear J;

% 150 pixels * 0.5mm = 7.5cm, 75mm corpus callosum is expected to be 5cm below the top of the skull and in the middle of the image
% we may change this when we start using acpcdetect
% when normcorrx2 runs
% the index cc(I, J) represents the filter's bottom-right point being ResampledAVW(I, J)
TopWMRow = min(I);
% expected centre of the CC is 
ExpectedCentreRow = round(TopWMRow + 75 ./ TemplatePixdims(2));
ExpectedCentreCol = round(size(ResampledAVW, 2) / 2);

ExpectedLowerRightI = round(ExpectedCentreRow + size(MatchedFilter, 1) / 2);
ExpectedLowerRightJ = round(ExpectedCentreCol + size(MatchedFilter, 2) / 2);
%ExpectedLowerRightI = round(ExpectedCentreRow - size(MatchedFilter, 1) / 2);
%ExpectedLowerRightJ = round(ExpectedCentreCol - size(MatchedFilter, 2) / 2);
%[ExpectedLowerRightI ExpectedLowerRightJ]
% T1
F = gaussian_filter_max_1d(2);
%keyboard;
WMSegSmoothed = imfilter(double(WMSeg), F(:), 'same', 'conv', 'replicate');
WMSegSmoothed = imfilter(WMSegSmoothed, F(:)', 'same', 'conv', 'replicate');

cc2 = normxcorr2(MatchedFilterProb, WMSegSmoothed);
% clf;
% 
% T = cc2(size(MatchedFilter, 1):end - size(MatchedFilter, 1) + 1, size(MatchedFilter, 2):end - size(MatchedFilter, 2) + 1);
% imshow(T, []);
% hold on;
% plot(ExpectedLowerRightJ, ExpectedLowerRightI, '*b');
%keyboard;
[cc2X, cc2Y] = meshgrid((((1:size(cc2, 2)) - ExpectedLowerRightJ)) .* TemplatePixdims(1), (((1:size(cc2, 1)) - ExpectedLowerRightI)) .* TemplatePixdims(2));

R = sqrt(cc2X .* cc2X + cc2Y .* cc2Y);
RWeighting = 0.5 * tanh(-(R - 50) ./ 10) + 0.5; 
%imshow(RWeighting, []);
%keyboard;
cc2R = cc2 .* RWeighting;
clear SEG;

cc2RRegMax = imregionalmax(cc2R, 4) & (cc2R > 0);
[cc2RRegMaxI, cc2RRegMaxJ] = find(cc2RRegMax);
[cc2RRegMaxIDX] = find(cc2RRegMax);
Rofcc2RRegMax = R(cc2RRegMax);
[~, I] = min(Rofcc2RRegMax);

yoffset2RReg = round(cc2RRegMaxI(I) - size(MatchedFilter, 1));
xoffset2RReg = round(cc2RRegMaxJ(I) - size(MatchedFilter, 2));

[~, SortedCCIDX] = sort(cc2R(:), 1, 'descend');
% only keep the regional maxima
SortedCCIDX = SortedCCIDX(ismember(SortedCCIDX, cc2RRegMaxIDX));
[SortedCCID, SortedCCJD] = ind2sub(size(cc2R), SortedCCIDX);

yoffset2R = round(SortedCCID - size(MatchedFilter, 1));
xoffset2R = round(SortedCCJD - size(MatchedFilter, 2));
% take out the nearest maxima to the centre, we are already doing this
I = yoffset2R == yoffset2RReg & xoffset2R == xoffset2RReg;
yoffset2R = yoffset2R(~I);
xoffset2R = xoffset2R(~I);
ExpectedCCRank = find(I);
%disp(['Expected location is of CC rank: ' num2str(ExpectedCCRank)]);
clear I;
%keyboard;

%FirstOffsetValid2R = find(yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1) & xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2), 1, 'first');
OffsetValid2R = find(...
	yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1) & ...
	xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2));

NumTries = 2;
Offsets = zeros(NumTries + 1, 2);
Offsets(1, :) = [yoffset2RReg, xoffset2RReg];
if(length(OffsetValid2R) < NumTries)
	disp('acpcdetect failed');
	ReturnCode = 1;
	return;
end
Offsets(2:NumTries + 1, :) = [yoffset2R(OffsetValid2R(1:NumTries)), xoffset2R(OffsetValid2R(1:NumTries))];

if(Offsets(1, 1) < 1 || Offsets(1, 1) + size(MatchedFilter, 1) - 1 > size(ResampledAVW, 1) || ...
   Offsets(1, 2) < 1 || Offsets(1, 2) + size(MatchedFilter, 2) - 1 > size(ResampledAVW, 2))
	Offsets = Offsets(2:NumTries + 1, :);
	disp('Centre cc invalid, removed');
end

%Offsets
% if we arent doing affine registration then just choose the start point
% that is closest to the centre which is the first offset pair
if(~DoLK)
	Offsets = Offsets(1, :);
end

InitialResampledAVW = cell(1, size(Offsets, 1));
MatchedFilterRemapped = cell(1, size(Offsets, 1));
LKParameters = cell(1, size(Offsets, 1));
LKCost = zeros(1, size(Offsets, 1));
TX = cell(1, size(Offsets, 1));
TY = cell(1, size(Offsets, 1));
InterpX = cell(1, size(Offsets, 1));
InterpY = cell(1, size(Offsets, 1));
TemplateLKIMG = cell(1, size(Offsets, 1));
TemplateProbLKIMG = cell(1, size(Offsets, 1));

% 
% T = cc2(size(MatchedFilter, 1):end - size(MatchedFilter, 1) + 1, size(MatchedFilter, 2):end - size(MatchedFilter, 2) + 1);
% imshow(T, []);
% hold on;
% plot(ExpectedLowerRightJ, ExpectedLowerRightI, '*b');
%G = cell(1, size(Offsets, 1));
FornixProbLKIMG = cell(1, size(Offsets, 1));
ResampledGroundCropped = cell(1, size(Offsets, 1));
ResampledAVWCropped = cell(1, size(Offsets, 1));
TemplateProbLKIMGCropped = cell(1, size(Offsets, 1));
FornixProbLKIMGCropped = cell(1, size(Offsets, 1));
OriginalOtsuMask = cell(1, size(Offsets, 1));
OtsuMask = cell(1, size(Offsets, 1));
TemplateOverlap = cell(1, size(Offsets, 1));
LKOtsuSEG = cell(1, size(Offsets, 1));

%DoLK
%keyboard;
for z = 1:size(Offsets, 1)
	[InitialResampledAVW{z}, ...
	MatchedFilterRemapped{z}, ...
	LKParameters{z}, LKCost(z), ...
	TX{z}, TY{z}, InterpX{z}, InterpY{z}, ...
	TemplateLKIMG{z}, TemplateProbLKIMG{z}, FornixProbLKIMG{z}, ...
	ResampledGroundCropped{z}, ResampledAVWCropped{z}, ...
	TemplateProbLKIMGCropped{z}, FornixProbLKIMGCropped{z}, ...
	OriginalOtsuMask{z}, OtsuMask{z}, ...
	TemplateOverlap{z}, LKOtsuSEG{z}] = do_lk_and_otsu(ResampledAVW, MatchedFilter, MatchedFilterProb, MatchedFilterFornixProb, ResampledGroundAVW, Offsets(z, :), DoLK);
	%TemplateOverlap{z}, LKOtsuSEG{z}] = do_lk_and_otsu(WMSegSmoothed, MatchedFilter, MatchedFilterProb, MatchedFilterFornixProb, ResampledGroundAVW, Offsets(z, :), DoLK);

end
% so sometimes the closest to the expected center fails and the centre needs to be the maximum correlation coefficient,
% so new strategy, try the closest and if this fails (OtsuMask is empty) then resort to the CC values
% find out which choice gives the best overlap on the template probabilities
%%
PCoverage = zeros(1, length(TemplateOverlap));

DiceCoefficients = zeros(1, length(TemplateOverlap));
Distances = cell(1, length(TemplateOverlap));
Angles = cell(1, length(TemplateOverlap));
SumDistances = zeros(1, length(TemplateOverlap));
SumAngles = zeros(1, length(TemplateOverlap));
TransformAngleX = zeros(1, length(TemplateOverlap));
TransformAngleY = zeros(1, length(TemplateOverlap));
LKCCSeg = cell(1, length(TemplateOverlap));

for z = 1:length(TemplateOverlap)
	PCoverage(z) = sum(TemplateProbLKIMGCropped{z}(:) .* OtsuMask{z}(:)) ./ sum(TemplateProbLKIMGCropped{z}(:));

	T = TemplateProbLKIMGCropped{z} > 0.3;
	LKCCSeg{z} = OtsuMask{z} & ~FornixProbLKIMGCropped{z} > 0.3;
	LKCCSeg{z} = bwareaopen(LKCCSeg{z}, 50);
	DiceCoefficients(z) = dices_coefficient(T, LKCCSeg{z});
	[Distances{z}, Angles{z}] = nearest_angle_distance(T, LKCCSeg{z});
	SumDistances(z) = sum(Distances{z});
	SumAngles(z) = sum(abs(Angles{z}));
%LKInitialOffsets
	TransformAngleX(z) = atan2(TY{z}(1, 2) - TY{z}(1, 1), TX{z}(1, 2) - TX{z}(1, 1));
	TransformAngleY(z) = atan2(TY{z}(2, 1) - TY{z}(1, 1), TX{z}(2, 1) - TX{z}(1, 1)) - pi / 2;
	clear T;
end
RigidFactor = abs(TransformAngleX - TransformAngleY);
ValidIDX = (TransformAngleX < pi / 3) & (TransformAngleY < pi / 3) & (RigidFactor < pi / 4);
if(ExpectedCCRank > 10)
	ValidIDX(1) = 0;
end
ValidIDXIDX = find(ValidIDX);

[MaxPCoverage, MaxPCoverageI] = max(PCoverage(ValidIDX));
[MaxDiceCoefficient, MaxDiceCoefficientI] = max(DiceCoefficients(ValidIDX));
% 	[MinHausdorffDistance, MinHausdorffDistanceI] = min(HausdorffDistances);
% 	[MinDeformationNorm, MinDeformationNormI] = min(DeformationNorm);
% 	[MinRotationAngle, MinRotationAngleI] = min(abs(RotationAngles));
[MinDistances, MinDistancesI] = min(SumDistances(ValidIDX));
[MinAngles, MinAnglesI] = min(SumAngles(ValidIDX));
%[MinTransformAngleX, MinTransformAngleXI] = min(TransformAngleX);
%[MinTransformAngleY, MinTransformAngleYI] = min(TransformAngleY);
[MinRigidFactor, MinRigidFactorI] = min(RigidFactor);
[MinLKCost, MinLKCostI] = min(LKCost(ValidIDX));
%disp(['Template overlap: ' num2str(MaxPCoverage)]);
%disp(['Best dice coefficient: ' num2str(MaxDiceCoefficient)]);
%disp(['Best Hausdorff distance: ' num2str(MinHausdorffDistance)]);

TestValues(1).text = 'PCoverage';	TestValues(1).bestidx = MaxPCoverageI;			TestValues(1).values = PCoverage;								TestValues(1).weight = 1;
TestValues(2).text = 'Dice';		TestValues(2).bestidx = MaxDiceCoefficientI;	TestValues(2).values = DiceCoefficients;						TestValues(2).weight = 1;
%TestValues(3).text = 'Hausdorff';	TestValues(3).bestidx = MinHausdorffDistanceI;	TestValues(3).values = HausdorffDistances;
TestValues(3).text = 'Angles';		TestValues(3).bestidx = MinAnglesI;				TestValues(3).values = SumAngles;								TestValues(3).weight = 1;
TestValues(4).text = 'Distances';	TestValues(4).bestidx = MinDistancesI;			TestValues(4).values = SumDistances;							TestValues(4).weight = 1;
TestValues(5).text = 'RigidFactor';	TestValues(5).bestidx = MinRigidFactorI;		TestValues(5).values = abs(TransformAngleX - TransformAngleY);	TestValues(5).weight = 1;
TestValues(6).text = 'LK Cost';		TestValues(6).bestidx = MinLKCostI;				TestValues(6).values = LKCost;									TestValues(6).weight = 3;

%BestIDXHist = hist([TestValues.bestidx], 1:length(TestValues));
%false
TestWeights = [TestValues.weight];

TestScores = zeros(1, length(ValidIDXIDX));
for z = 1:length(ValidIDXIDX)
	TestScores(z) = sum(TestWeights([TestValues.bestidx] == z));
end
%TestScores
%ValidIDXIDX
%keyboard;
if(~any(ValidIDX))
	BestAlignmentIDX = 1;
else
	[~, BestAlignmentIDX] = max(TestScores);
	BestAlignmentIDXOriginal = ValidIDXIDX(BestAlignmentIDX);
end

DoCCGraphics = true;
if(DoCCGraphics)
	clf;
	SR = 2;
	SC = size(Offsets, 1);
	%keyboard;
	%clf;
	iptsetpref('ImshowAxesVisible', 'off');
	AX = zeros(SR, SC);
	for z = 1:size(Offsets, 1)

	%	title({'Normxcorr result, red based on maximum cc,', 'blue based on minimum distance from expected'});
		[~, LOC] = ismember(z, ValidIDXIDX);

		AX(1, z) = subplot(SR, SC, z);
		imshow(ResampledAVW, []);
		hold on;
		SZ = size(MatchedFilter);

		rectangle('Position', [Offsets(z, 2), Offsets(z, 1), SZ(2), SZ(1)], 'EdgeColor', 'g', 'LineWidth', 2);
% 			T = line([Offsets(z, 2), Offsets(z, 2) + SZ(2)], [Offsets(z, 1), Offsets(z, 1)]); set(T, 'Color', 'b');
% 			T = line([Offsets(z, 2), Offsets(z, 2) + SZ(2)], [Offsets(z, 1) + SZ(1), Offsets(z, 1) + SZ(1)]); set(T, 'Color', 'b');
% 			T = line([Offsets(z, 2), Offsets(z, 2)], [Offsets(z, 1), Offsets(z, 1) + SZ(1)]); set(T, 'Color', 'b');
% 			T = line([Offsets(z, 2) + SZ(2), Offsets(z, 2) + SZ(2)], [Offsets(z, 1), Offsets(z, 1) + SZ(1)]); set(T, 'Color', 'b');

		LineProps = {'Color', 'r', 'LineWidth', 2};
		T = line([TX{z}(1, 1) TX{z}(end, 1)], [TY{z}(1, 1), TY{z}(end, 1)]); set(T, LineProps{:});
		T = line([TX{z}(1, end) TX{z}(end, end)], [TY{z}(1, end), TY{z}(end, end)]); set(T, LineProps{:});
		T = line([TX{z}(1, 1) TX{z}(1, end)], [TY{z}(1, 1), TY{z}(1, end)]); set(T, LineProps{:});
		T = line([TX{z}(end, 1) TX{z}(end, end)], [TY{z}(end, 1), TY{z}(end, end)]); set(T, LineProps{:});
		title(roman_label(z), 'FontSize', 20);
		%[MaxPCoverage, MaxPCoverageI] = max(PCoverage);
		%[MaxDiceCoefficient, MaxDiceCoefficientI] = max(DiceCoefficients);
		%[MinHausdorffDistance, MinHausdorffDistanceI] = min(HausdorffDistances);

		AX(2, z) = subplot(SR, SC, z + size(Offsets, 1));
		T = zeros(size(LKCCSeg{z}));
		S = TemplateProbLKIMGCropped{z} > 0.3;
		T(~S & LKCCSeg{z}) = 2;
		T(S & ~LKCCSeg{z}) = 3;
		T(S & LKCCSeg{z}) = 4;
		imshow(T, []);
		clear S T;

		C = cell(1, length(TestValues) + 1);
		for TestIDX = 1:length(TestValues)
			if(TestValues(TestIDX).values(z) == floor(TestValues(TestIDX).values(z)))
				FormatString = '%d';
			else
				FormatString = '%.2f';
			end
			C{TestIDX} = ['{\it' TestValues(TestIDX).text '}: ' num2str(TestValues(TestIDX).values(z), FormatString)];
			if(TestValues(TestIDX).bestidx == LOC)
				C{TestIDX} = ['{\bf' C{TestIDX} '}'];
			end
		end
		if(ValidIDX(z) == false)
			TextColor = [0.5, 0.5, 0.5];
			C{end} = 'Invalid';
		elseif BestAlignmentIDX == LOC
			TextColor = 'k';
			C = C(1:end - 1);
		else
			TextColor = [0.5, 0.5, 0.5];
			C = C(1:end - 1);
		end

		text(0, -0.01, C, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Color', TextColor);
		clear C;
% 		
% 			subplot(SR, SC, z + size(Offsets, 1));
% 			%imshow(OriginalOtsuMask{z}, []);
% 			imshow(LKOtsuSEG{z}, []);
% 			hold on;
% 			[~, CC] = contour(TemplateProbLKIMGCropped{z}, [0.3, 0.3]);
% 			%[~, CC] = contour(NewPSI, [0.3, 0.3]);
% 			set(CC, 'Color', 'r');
% 			xlabel('Otsu mask, template is red contour');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
		AXPos = get(AX(2, z), 'Position');
		AXPos(2) = AXPos(2) + 0.125;
		set(AX(2, z), 'Position', AXPos);
	end
	LegHandle = text(1.01, 0.5, {'   Dark Grey - {\itLKCCSeg}', '   Light Grey - {\itPCCMask}', '   White - {\itLKCCSeg} \cap {\itPCCMask}'}, ...
		'Units', 'Normalized', ...
		'Color', 'k', ...
		'HorizontalAlignment', 'left', ...
		'VerticalAlignment', 'middle', ...
		'FontSize', 14, ...
		'Parent', AX(2, size(Offsets, 1)));

	LegHandle = text(1.01, 0.5, {'Affine Alignment', '  Green - Initial', '  Red - After LKT'}, ...
		'Units', 'Normalized', ...
		'Color', 'k', ...
		'HorizontalAlignment', 'left', ...
		'VerticalAlignment', 'middle', ...
		'FontSize', 14, ...
		'Parent', AX(1, size(Offsets, 1)));

	for z = 1:size(AX, 2)
		N = 0.05;
		AXPos = get(AX(1, z), 'Position');
		AXPos(1) = AXPos(1) - N * z;
		set(AX(1, z), 'Position', AXPos);
		AXPos = get(AX(2, z), 'Position');
		AXPos(1) = AXPos(1) - N * z;
		set(AX(2, z), 'Position', AXPos);

	end

%	if MainFittingIteration == 1
		OutputFile = fullfile(OutputDir, 'cc', [OutputPrefix '_cc_choice.png']);
%	else
%		OutputFile = fullfile(OutputDir, 'cc', [OutputPrefix '_cc_fallback.png']);
%	end
	OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-002-LKCC.eps']);
	%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
	
	keyboard;
	cc_seg_save_figure_paper_eps(OutputFile);
	cc_seg_save_figure_paper_png(OutputFile);
% 		FigPos = fullscreen_fig_pos;
% 		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% 		exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 		IMG = imread(OutputFile);
% 		IMG = imautocropwhite(IMG);
% 		imwrite(IMG, OutputFile);
end

if(~all(isinf(SumDistances)) && any(ValidIDX))
	InitialResampledAVW = InitialResampledAVW(ValidIDX);
	MatchedFilterRemapped = MatchedFilterRemapped(ValidIDX);
	LKParameters = LKParameters(ValidIDX);
	TX = TX(ValidIDX);
	TY = TY(ValidIDX);
	InterpX = InterpX(ValidIDX);
	InterpY = InterpY(ValidIDX);
	TemplateLKIMG = TemplateLKIMG(ValidIDX);
	TemplateProbLKIMG = TemplateProbLKIMG(ValidIDX);
	FornixProbLKIMG = FornixProbLKIMG(ValidIDX);
	ResampledAVWCropped = ResampledAVWCropped(ValidIDX);
	ResampledGroundCropped = ResampledGroundCropped(ValidIDX);
	TemplateProbLKIMGCropped = TemplateProbLKIMGCropped(ValidIDX);
	FornixProbLKIMGCropped = FornixProbLKIMGCropped(ValidIDX);
	OriginalOtsuMask = OriginalOtsuMask(ValidIDX);
	OtsuMask = OtsuMask(ValidIDX);
	TemplateOverlap = TemplateOverlap(ValidIDX);
	LKCCSeg = LKCCSeg(ValidIDX);

	BestLKI = BestAlignmentIDX;
	%keyboard;
	InitialResampledAVW = InitialResampledAVW{BestLKI};
	MatchedFilterRemapped = MatchedFilterRemapped{BestLKI};
	LKParameters = LKParameters{BestLKI};
	TX = TX{BestLKI};
	TY = TY{BestLKI};
	InterpX = InterpX{BestLKI};
	InterpY = InterpY{BestLKI};
	TemplateLKIMG = TemplateLKIMG{BestLKI};
	TemplateProbLKIMG = TemplateProbLKIMG{BestLKI};
	FornixProbLKIMG = FornixProbLKIMG{BestLKI};
	ResampledAVWCropped = ResampledAVWCropped{BestLKI};
	ResampledGroundCropped = ResampledGroundCropped{BestLKI};
	TemplateProbLKIMGCropped = TemplateProbLKIMGCropped{BestLKI};
	FornixProbLKIMGCropped = FornixProbLKIMGCropped{BestLKI};
	OriginalOtsuMask = OriginalOtsuMask{BestLKI};
	OtsuMask = OtsuMask{BestLKI};
	LKCCSeg = LKCCSeg{BestLKI};
	TemplateOverlap = TemplateOverlap{BestLKI};
else
%	if MainFittingIteration == 2
		disp('Something went wrong, check your file');
		ReturnCode = 1;
		return;
%	end
end

LKCCSeg = LKCCSeg & (TemplateProbLKIMGCropped > 0.5);

MeanIntensityInCCSeg = mean(ResampledAVWCropped(LKCCSeg));
VarIntensityInCCSeg = var(ResampledAVWCropped(LKCCSeg));

GaussianProb = normpdf(ResampledAVWCropped, MeanIntensityInCCSeg, sqrt(VarIntensityInCCSeg));
GaussianProb = GaussianProb ./ max(GaussianProb(:));
GaussianProb(ResampledAVWCropped > MeanIntensityInCCSeg) = 1;

SIGMA = 0.25;

[GaussianFilter, GaussianFilterDeriv] = gaussian_filter_max_1d(SIGMA);

aSmooth=imfilter(GaussianProb, GaussianFilter, 'conv', 'replicate', 'same');   % run the filter across rows
aSmooth=imfilter(aSmooth, GaussianFilter','conv','replicate', 'same'); % and then across columns

%apply directional derivatives
ax = imfilter(aSmooth, GaussianFilterDeriv(:)', 'conv','replicate');
ay = imfilter(aSmooth, GaussianFilterDeriv(:), 'conv','replicate');

GaussianGradMAG = sqrt(ax .* ax + ay .* ay);
clear ax ay aSmooth;
GaussianGradMAG = GaussianGradMAG ./ max(GaussianGradMAG(:));

PenaltyA = min(FornixProbLKIMGCropped * 50, 1) .* GaussianProb;
PenaltyB = GaussianGradMAG;
PenaltyC = exp(-2 .* GaussianProb .* (TemplateProbLKIMGCropped + 2));
PenaltyA = PenaltyA ./ max(PenaltyA(:));
PenaltyB = PenaltyB ./ max(PenaltyB(:));
PenaltyC = PenaltyC ./ max(PenaltyC(:));

PenaltyImage = PenaltyA + 4 * PenaltyB + 2 * PenaltyC;

[StrongI, StrongJ] = find(LKCCSeg);
CCSeg = bwselect(PenaltyImage < 0.8, StrongJ, StrongI, 8);
clear StrongI StrongJ;

CCSeg = bwareaopen(CCSeg, 200);

CCSeg = bwfill(CCSeg, 'holes', 8);
CCSegCC = bwconncomp(CCSeg, 8);
NumRegionsCCSeg = CCSegCC.NumObjects;
if(NumRegionsCCSeg > 1)

	T = bw_join_mst_by_boundaries(CCSeg);
	CCSeg = CCSeg | (T & GaussianProb > 0.3);
	clear T;
end

CCSeg = bwfill(CCSeg, 'holes', 8);
CCSegCC = bwconncomp(CCSeg, 8);
NumRegionsCCSeg = CCSegCC.NumObjects;

if(NumRegionsCCSeg > 1)
	try
		NewCCSegJoiningSegments = bw_join_mst_by_boundaries(CCSeg);
	catch e
		disp('acpcdetect failed');
		ReturnCode = 1;
		return;
	end
	CCSeg = CCSeg | NewCCSegJoiningSegments;
end
	
%keyboard;
clear Overlaps DMaskL L;

SE = strel('disk', 4);
%DMask = imclose(DMask, SE);
CCSeg = imclose(CCSeg, SE);

CCSegCC = bwconncomp(CCSeg);
CCSegL = labelmatrix(CCSegCC);
CCSegAreas = cellfun('length', CCSegCC.PixelIdxList);
[~, L] = max(CCSegAreas);
if(isempty(L))
	ReturnCode = 1;
	return;
end
CCSeg = (CCSegL == L);

% DMask = bwfill(DMask, 'holes');
CCSeg = bwfill(CCSeg, 'holes');
BW = PenaltyImage < 1 | CCSeg;
I = find(CCSeg(:), 1, 'first');
[ID, JD] = ind2sub(size(CCSeg), I);
BWS = bwselect(BW, JD, ID, 8);
AddedByReconstruct = BWS & ~CCSeg;
clear I ID JD BWS BW;
CCSeg = CCSeg | AddedByReconstruct;

CCSeg = imclose(CCSeg, strel('disk', 3));

FinalSeg = CCSeg;

T = ResampledAVWCropped .* FinalSeg;
T(~FinalSeg) = -Inf;
SE = strel('disk', 3);
T = imdilate(T, SE);
T(~FinalSeg) = Inf;
TopHat = T - ResampledAVWCropped;
TopHat(~FinalSeg) = 0;

T = ResampledAVWCropped .* FinalSeg;
T(~FinalSeg) = Inf;
SE = strel('disk', 3);
T = imerode(T, SE);
T(~FinalSeg) = Inf;
BottomHat = ResampledAVWCropped - T;
BottomHat(~FinalSeg) = 0;

MorphGrad = BottomHat + TopHat;

OutlierPixelsStrong = MorphGrad > (median(MorphGrad(FinalSeg)) * 3);
OutlierPixelsWeak = MorphGrad > (median(MorphGrad(FinalSeg)) * 2);
clear MorphGrad;

if(any(OutlierPixelsStrong(:)))
	[OutlierPixelsStrongR, OutlierPixelsStrongC] = find(OutlierPixelsStrong);
	OutlierPixels = bwselect(OutlierPixelsWeak, OutlierPixelsStrongC, OutlierPixelsStrongR, 8);
else
	OutlierPixels = false(size(OutlierPixelsStrong));
end

R = regionprops(double(FinalSeg), 'BoundingBox');
LeftColumn = floor(R.BoundingBox(1) + R.BoundingBox(3) * 2 / 3);
OutlierPixels(:, 1:LeftColumn) = 0;
clear R LeftColumn;

% dilate the outlier pixels
if(~any(OutlierPixels(:)))
	FinalSegOpenedM	= FinalSeg;
else
	SE = strel('disk', 3);
	
	OutlierPixelsDilated = imdilate(OutlierPixels, SE);
	FinalSegWithout = FinalSeg & ~OutlierPixelsDilated;
	T = imdilate(OutlierPixelsDilated, ones(3, 3));
	FinalSegWithoutOutlierBorders = T & FinalSegWithout;
	FinalSegOpened = imdilate(FinalSegWithoutOutlierBorders, SE) | FinalSegWithout;%imopen(FinalSegFirstCull, strel('disk', 3));
    clear FinalSegWithoutOutlierBorders;
	FinalSegOpened = imclose(FinalSegOpened, strel('disk', 2));
	FinalSegOpened(~FinalSeg) = 0;
	FinalSegOpenedCC = bwconncomp(FinalSegOpened);
	FinalSegOpenedAreas = cellfun('length', FinalSegOpenedCC.PixelIdxList);
	FinalSegOpenedL = labelmatrix(FinalSegOpenedCC);
	if(sum(FinalSegOpenedAreas > 200) > 1)
		[JoiningSegments] = bw_join_mst_by_boundaries(FinalSegOpened);
		FinalSegOpened = FinalSegOpened | JoiningSegments;
		FinalSegOpenedCC = bwconncomp(FinalSegOpened);
		FinalSegOpenedAreas = cellfun('length', FinalSegOpenedCC.PixelIdxList);
		FinalSegOpenedL = labelmatrix(FinalSegOpenedCC);
	end
	
	[~, I] = max(FinalSegOpenedAreas);
	FinalSegOpenedM = FinalSegOpenedL == I;
	FinalSegOpenedM = imclose(FinalSegOpenedM, SE);
	FinalSegOpenedM = bwfill(FinalSegOpenedM, 'holes');
	
end

FinalSegOpenedMOpened = imopen(FinalSegOpenedM, strel('disk', 1));
FinalSegOpenedMOpenedLost = ~FinalSegOpenedMOpened & FinalSegOpenedM;

FinalSegOpenedMOpenedLostCC = bwconncomp(FinalSegOpenedMOpenedLost, 8);
FinalSegOpenedMThickened = FinalSegOpenedMOpened;
for z = 1:length(FinalSegOpenedMOpenedLostCC.PixelIdxList)
	T = FinalSegOpenedM;
	T(FinalSegOpenedMOpenedLostCC.PixelIdxList{z}) = 0;
	TCC = bwconncomp(T);
	if(TCC.NumObjects > 1)
		T = false(size(FinalSegOpenedMOpened));
		T(FinalSegOpenedMOpenedLostCC.PixelIdxList{z}) = 1;
		FinalSegOpenedMThickened = FinalSegOpenedMThickened | imdilate(T, strel('square', 3));
	end
	clear T TCC;
end

FinalSegOpenedMThickened = imdilate(FinalSegOpenedM, strel('disk', 1));
FinalSegOpenedMThickened = bwfill(FinalSegOpenedMThickened, 'holes', 8);

FinalSeg = imdilate(FinalSeg, strel('disk', 1));
FinalSeg = bwfill(FinalSeg, 'holes', 8);

OutputS.InitialSeg = CCSeg;
OutputS.FinalSeg = FinalSeg;
OutputS.FinalSegArtefactsRemoved = FinalSegOpenedMThickened;
OutputS.IMG = ResampledAVWCropped;
OutputS.TemplatePixdims = TemplatePixdims;

if(GroundTruthGiven)
	OutputS.GroundSeg = ResampledGroundCropped > 0;
	OutputS.GroundSeg = bwfill(OutputS.GroundSeg, 'holes', 8);
	OutputS.InitialDice = dices_coefficient(OutputS.InitialSeg, OutputS.GroundSeg);
	OutputS.FinalDice = dices_coefficient(OutputS.FinalSegArtefactsRemoved, OutputS.GroundSeg);
	
	FinalSegBoundaries = bwboundaries(OutputS.FinalSeg, 8);
	GroundSegBoundaries = bwboundaries(OutputS.GroundSeg, 8);
	if(length(FinalSegBoundaries) > 1 || length(GroundSegBoundaries) > 1)
		error('There should only be one boundary in the final and ground segmentations');
	end
	
	OutputS.InitialHaussdorf = hausdorff_distance(OutputS.InitialSeg, OutputS.GroundSeg, TemplatePixdims);
	OutputS.FinalHaussdorf = hausdorff_distance(OutputS.FinalSegArtefactsRemoved, OutputS.GroundSeg, TemplatePixdims);
end

save(SegMatFileName, '-struct', 'OutputS');
ReturnCode = 0;

function [H] = hausdorff_distance(A, B, Pixdims)

BoundaryA = bwboundaries(A);
BoundaryB = bwboundaries(B);

if(~isempty(BoundaryA) && ~isempty(BoundaryB))
	BoundaryA = cat(1, BoundaryA{:});
	BoundaryB = cat(1, BoundaryB{:});
	BoundaryA = bsxfun(@times, BoundaryA, Pixdims(:)');
	BoundaryB = bsxfun(@times, BoundaryB, Pixdims(:)');
	DX = bsxfun(@minus, BoundaryA(:, 2), BoundaryB(:, 2)');
	DY = bsxfun(@minus, BoundaryA(:, 1), BoundaryB(:, 1)');
	D = sqrt(DX .* DX + DY .* DY);
	H = max(max(min(D)), max(min(D, [], 2)));
else
	H = Inf;
end

function [Distances, Angles, SmoothedBoundaryA, SmoothedBoundaryB] = nearest_angle_distance(A, B)
% computes the sum of the differences in boundary vectors
% so for each boundary point in A, computes the angle between the tangent
% vector and the tangent vector of the nearest boundary point in B

BoundaryA = bwboundaries(A);
BoundaryB = bwboundaries(B);

if(~isempty(BoundaryA) && ~isempty(BoundaryB))
	% smooth the boundaries first before computing tangents
	SmoothedBoundaryA = BoundaryA;
	SmoothedBoundaryB = BoundaryB;
	
	for z = 1:length(SmoothedBoundaryA)
		SmoothedBoundaryA{z} = imfilter(SmoothedBoundaryA{z}, [1, 1, 1]' / 3, 'same', 'conv', 'circular');
	end
	for z = 1:length(SmoothedBoundaryB)
		SmoothedBoundaryB{z} = imfilter(SmoothedBoundaryB{z}, [1, 1, 1]' / 3, 'same', 'conv', 'circular');
	end
	
	TangentA = SmoothedBoundaryA;
	TangentB = SmoothedBoundaryB;
	
	for z = 1:length(TangentA)
		TangentA{z} = imfilter(TangentA{z}, [1, 0, -1]', 'same', 'conv', 'circular');
	end
	for z = 1:length(TangentB)
		TangentB{z} = imfilter(TangentB{z}, [1, 0, -1]', 'same', 'conv', 'circular');
	end
	
	SmoothedBoundaryA = cat(1, SmoothedBoundaryA{:});
	SmoothedBoundaryB = cat(1, SmoothedBoundaryB{:});
	TangentA = cat(1, TangentA{:});
	TangentB = cat(1, TangentB{:});
	
	TangentANorm = sqrt(sum(TangentA .* TangentA, 2));
	TangentBNorm = sqrt(sum(TangentB .* TangentB, 2));
	
	TangentA = bsxfun(@rdivide, TangentA, TangentANorm);
	TangentB = bsxfun(@rdivide, TangentB, TangentBNorm);
	
	DX = bsxfun(@minus, SmoothedBoundaryA(:, 2), SmoothedBoundaryB(:, 2)');
	DY = bsxfun(@minus, SmoothedBoundaryA(:, 1), SmoothedBoundaryB(:, 1)');
	D = sqrt(DX .* DX + DY .* DY);
	[Distances, DistanceIDX] = min(D, [], 2);
	
	% perform dot product
	T = sum(TangentA .* TangentB(DistanceIDX, :), 2);
	% clip
	T = min(T, 1);
	T = max(T, -1);
	Angles = acos(T);
else
	Distances = Inf;
	Angles = Inf;
end
%keyboard;

function [D] = dices_coefficient(A, B)

D = 2 .* sum(A(:) & B(:)) ./ (sum(A(:)) + sum(B(:)));

function [JoiningSegments] = bw_join_mst_by_boundaries(BW)

[B] = bwboundaries(BW, 8);
if(length(B) == 1)
	JoiningSegments = false(size(BW));
else

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
		% this is because we calculate the distances based on increasing indices of edges
		% if the edges are returned with a greater index on the left it
		% will create errors downstream, so sort them
		edges = sort(edges, 2); 
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
		I = sub2ind(size(BW), IY, IX);
		I = unique(I);
		T = false(size(BW));
		T(I) = 1;
		Angle = atan2(EndY - StartY, EndX - StartX);
		Angle = mod(Angle + 2 * pi, 2 * pi);

		R = [5, 5] + 7 * abs([cos(Angle), sin(Angle)]);

		AngleWeighting = -0.9 * cos(2 * (Angle + 45 * pi / 180));

		SQ = sqrt(R(1) * R(2)) * AngleWeighting;
		SIGMA = [R(1), SQ; SQ, R(2)];
		[F, HalfMaximum] = gaussian_fwhm2d(SIGMA);

		JoiningSegments = imdilate(T, double(F > HalfMaximum));

		clear I IX IY n ArcLength EndY StartY EndX StartX F T HalfMaximum SQ Angle AngleWeighting R;
	end
end

function [AVWMidSag] = extract_mid_slice(AVW)

if(mod(size(AVW, 2), 2) == 0)
	XI = size(AVW, 2) / 2;
	XFrac = 0.5;
	AVWMidSag = AVW(size(AVW, 1):-1:1, XI, size(AVW, 3):-1:1) .* (1 - XFrac) + ...
	AVW(size(AVW, 1):-1:1, XI + 1, size(AVW, 3):-1:1) .* XFrac;
	clear XI XFrac MidsagSliceIDX;
else
	XI = ceil(size(AVW, 2) / 2);
	AVWMidSag = AVW(size(AVW, 1):-1:1, XI, size(AVW, 3):-1:1);
	clear XI;
end
AVWMidSag = squeeze(AVWMidSag)';
AVWMidSag = double(AVWMidSag);

function [InitialResampledAVW, ...
	MatchedFilterRemapped, ...
	LKParameters, LKCost, ...
	TX, TY, InterpX, InterpY, ...
	TemplateLKIMG, TemplateProbLKIMG, FornixProbLKIMG, ...
	ResampledGroundCropped, ResampledAVWCropped, ...
	TemplateProbLKIMGCropped, FornixProbLKIMGCropped, ...
	OriginalOtsuMask, OtsuMask, ...
	TemplateOverlap, OtsuSEG] = do_lk_and_otsu(ResampledAVW, MatchedFilter, MatchedFilterProb, MatchedFilterFornixProb, ResampledGroundAVW, real_total_offset, DoLK)

if nargin < 7
	DoLK = true;
end
%real_total_offset = [118 116];
InitialResampledAVW = ResampledAVW(real_total_offset(1):real_total_offset(1) + size(MatchedFilter, 1) - 1, real_total_offset(2):real_total_offset(2) + size(MatchedFilter, 2) - 1);

MatchedFilterRemapped = immatchhist(MatchedFilter, InitialResampledAVW);
%keyboard;
if(DoLK)
	[LKParameters, LKCost] = lk_weighted_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, MatchedFilter ./ max(MatchedFilter(:)), [0 0 0 0 real_total_offset(2) real_total_offset(1)], 100);
else
	[LKParameters, LKCost] = lk_weighted_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, MatchedFilter ./ max(MatchedFilter(:)), [0 0 0 0 real_total_offset(2) real_total_offset(1)], 0);
	%LKParameters = [0 0 0 0 real_total_offset(2) real_total_offset(1)]';
	%LKCost = 0;
end
[TX, TY, InterpX, InterpY] = coords_template_lk_img(LKParameters, ResampledAVW, MatchedFilter);
TemplateLKIMG = interp2(MatchedFilter, InterpX, InterpY, 'linear', 0);
TemplateProbLKIMG = interp2(MatchedFilterProb, InterpX, InterpY, 'linear', 0);
FornixProbLKIMG = interp2(MatchedFilterFornixProb, InterpX, InterpY, 'linear', 0);

%% find the bounding box of the template warped to the image
% from regionprops
[listI, listJ] = find(TemplateLKIMG > 0);
list = [listJ, listI];
min_corner = min(list,[],1) - 0.5;
max_corner = max(list,[],1) + 0.5;
TemplateLKIMGBoundingBox = [min_corner (max_corner - min_corner)];
clear min_corner max_corner list listI listJ;

ResampledAVWCropped = imcrop(ResampledAVW, TemplateLKIMGBoundingBox);

if(~isempty(ResampledGroundAVW))
	ResampledGroundCropped = imcrop(ResampledGroundAVW, TemplateLKIMGBoundingBox);
else
	ResampledGroundCropped = [];
end
TemplateProbLKIMGCropped = imcrop(TemplateProbLKIMG, TemplateLKIMGBoundingBox);
FornixProbLKIMGCropped = imcrop(FornixProbLKIMG, TemplateLKIMGBoundingBox);

[THRESH, OtsuSEG] = robust_otsu2(ResampledAVWCropped, [0.05, 0.98]);
OtsuSEG(ResampledAVWCropped > THRESH(end)) = 3;
OtsuMask = (OtsuSEG == 3);

OriginalOtsuMask = OtsuMask;
FornixMaskIdx = find((FornixProbLKIMGCropped * 50) > 0.5);

OtsuMaskCC = bwconncomp(OtsuMask);
OtsuMaskL = labelmatrix(OtsuMaskCC);
OtsuMaskR = regionprops(OtsuMaskCC, 'Area', 'PixelIdxList');

TemplateOverlap = zeros(length(OtsuMaskR), 1);
% keep all regions that have high overlap with template or big area
for z = 1:length(OtsuMaskR)
	IDX = setdiff(OtsuMaskR(z).PixelIdxList, FornixMaskIdx);
	TemplateOverlap(z) = sum(TemplateProbLKIMGCropped(IDX)) ./ length(IDX);
	clear IDX;
end
clear FornixMaskIdx;
I = find(TemplateOverlap > 0.01);
if(~isempty(I))
	OtsuMask = ismember(OtsuMaskL, I);
	OtsuMask = bwareaopen(OtsuMask, 200);
	clear Junk OtsuMaskL OtsuMaskR I;
else
	disp('Current XCORR and LK failed');
	OtsuMask = false(size(OtsuMaskL));
end

function [M] = mean_non_nan(A)

M = zeros(1, size(A, 2));
NumNonNaN = sum(~isnan(A));
I = NumNonNaN > 0;

if(any(I))
	T = A(:, I);
	T(isnan(T)) = 0;

	M(I) = sum(T) ./ NumNonNaN(I);
end
