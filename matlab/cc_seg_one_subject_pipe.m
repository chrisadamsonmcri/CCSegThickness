function [varargout] = cc_seg_one_subject_pipe(FileNameOrStruct, GroundTruthFile, OutputDir, OutputPrefix)

% cc_seg_one_subject_pipe(NIIFileName)
% cc_seg_one_subject_pipe(Data)
%
% DESCRIPTION
%	Performs corpus callosum segmentation the using illusory contours level set method and various mathematical morphology operations.
%
% PARAMETERS
%	First argument
%	NIIFileName (char): expected to be a NIFTI file
%	Data (struct) with fields:
%		.AVW: midsagittal slice
%		.PixDims: pixel dimensions (X, Y) of midasgittal slice

% checks to see if we have a char array (so a nifti file) or a struct (CC image is a matlab array obtained in another function)
% what needs to come out of this
% AVWMidSag is the 2D midsagittal volume
% NIIPixdims [2]: X and Y dimensions of the image AVW

if(ischar(FileNameOrStruct))
	disp('We are processing a file');
	InFile = FileNameOrStruct;
	[NII, AVW] = load_nii(InFile);
	if(ismatrix(AVW))
		disp('A 2D file');
		% 2D image
		AVWMidSag = double(AVW);
		T = NII.hdr.dime.pixdim(2:4);
		NIIPixdims = T(NII.hdr.dime.dim(2:4) > 1);
		clear T;
	elseif(ndims(AVW) == 3)
		disp('A 3D file');
		[pathstr, name, ext] = fileparts(FileNameOrStruct);
		OutputMAT = fullfile(pathstr, [name '_art.mat']);
		clear ext;
		clear pathstr name;
		
		if(exist(OutputMAT, 'file') == 2)
			disp('Art already run, loading from file');
			load(OutputMAT);
		else
			disp('Art running');
			% 3D image, pipe through acpcdetect and extract middle slice
			ARTTempDirectory = tempname;
			if(exist(ARTTempDirectory, 'file') == 2)
				delete(ARTTempDirectory);
			elseif(exist(ARTTempDirectory, 'dir') == 7)
				rmdir(ARTTempDirectory);
			end

			[~, ~, ~] = mkdir(ARTTempDirectory);

			% use fslswapdim to put the file into canonical orientation
			[~, result] = system(['export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; fslorient -getorient ' InFile]);

			%disp(result);

			switch(deblank(result))
				case 'NEUROLOGICAL'
					OrientString = 'LR PA IS';
				case 'RADIOLOGICAL'
					OrientString = 'RL PA IS';
			end

			NIIFileForART = fullfile(ARTTempDirectory, 'in.nii');

			[~, ~] = system(['export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; fslswapdim ' InFile ' ' OrientString ' ' NIIFileForART]);

			[NII, AVW] = load_nii(NIIFileForART);

			% run cropping now rather than in python script
			
			T = double(sort(AVW(:)));
			T = T(floor(length(T) * 0.02):ceil(length(T) * 0.98));
			
			OriginalRange = max(T(:)) - min(T(:));
			OriginalMin = min(T(:));
			
			T = (T - OriginalMin) ./ OriginalRange;
			T = uint8(round(T * 255));
			
			OtsuThresh = otsu1(T);
			clear T;
			OtsuSegOfAVW = false(size(AVW));
			OtsuSegOfAVW(AVW > OtsuThresh / 255 * OriginalRange + OriginalMin) = 1;
			
			OtsuSegOfAVWCC = bwconncomp(OtsuSegOfAVW, 26);
			OtsuSegOfAVWAreas = cellfun('length', OtsuSegOfAVWCC.PixelIdxList);
			[~, I] = max(OtsuSegOfAVWAreas);
			
			[ID, JD, KD] = ind2sub(size(AVW), OtsuSegOfAVWCC.PixelIdxList{I});
			BoundingBoxI = [min(ID), max(ID)];
			BoundingBoxJ = [min(JD), max(JD)];
			BoundingBoxK = [min(KD), max(KD)];
			
			BoundingBoxK(1) = floor(BoundingBoxK(2) - 180 / NII.hdr.dime.pixdim(4));
			BoundingBoxK(1) = max(BoundingBoxK(1), 1);
			%SliceDownTwentyCM = int(N.floor(MaxSlice - 180 / InFileNII.header['pixdim'][3]))
			%if SliceDownTwentyCM < 0:
				%SliceDownTwentyCM = 0 
			AVW = AVW(BoundingBoxI(1):BoundingBoxI(2), BoundingBoxJ(1):BoundingBoxJ(2), BoundingBoxK(1):BoundingBoxK(2));
			clear BoundingBoxI BoundingBoxJ BoundingBoxK ID JD KD I OtsuSegOfAVWAreas OtsuSegOfAVWCC;
			%I = find(OtsuSegOfAVW);
			%[ID, JD, KD] = 
			%BoundingBoxK(2) = [Boundimin(KD), max(KD)];
			%keyboard;
						
			% scale between 0 and 1000
			AVW = double(AVW);
			M = max(AVW(:));
			if(M > 32767)
				AVW = (AVW - min(AVW(:))) ./ (max(AVW(:)) - min(AVW(:)));
				AVW = int16(round(AVW * 1000));
			else
				AVW = int16(AVW);
			end

			NII.img = permute(flipdim(AVW, 1), [2 1 3]);
			save_nii(NII, NIIFileForART);

			%keyboard;
			% make sure the file is in the correct data type for ART

			%[~, ~] = system(['export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; fslmaths ' NIIFileForART ' ' NIIFileForART ' -odt short']);
			%keyboard;

			if(exist([NIIFileForART '.gz'], 'file') == 2)
				%gunzip([NIIFileForART '.gz']);
				delete([NIIFileForART '.gz']);
			end

			% now perform cropping

% 			CurrentDir = pwd;
% 			cd(ARTTempDirectory);
% 
% 			if(exist('/home/addo/Dropbox/matlab/cc_seg/reorient_crop_one_subject.sh', 'file') == 2)
% 				ScriptName = '/home/addo/Dropbox/matlab/cc_seg/reorient_crop_one_subject.sh';
% 			else
% 				ScriptName = '/home/addo/matlab/cc_seg/reorient_crop_one_subject.sh';
% 			end
% 
% 			CommandString = ['export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; ' ScriptName ' in'];
% 			disp(CommandString);
% 			[status, result] = system(CommandString);
% 			[NIIOtsu, AVWOtsu] = load_nii('in_otsu');
% 			[NIIFinal, AVWFinal] = load_nii('in_maxlabel_cropped');
% 			keyboard;
% 			
% 			%%
% 			C = 71;
% 			subplot 121;
% 			imshow(OtsuSegOfAVW(:, :, C), [])
% 			subplot 122;
% 			imshow(AVWOtsu(:, :, C), [])
% 			
% 			%%
% 			cd(CurrentDir);
% 			keyboard;
	% 		[match] = regexp(InFile, '.*\.nii\.gz$', 'match');
	% 		if(~isempty(match))
	% 			CompressedInput = true;
	% 			NIIFileForARTPrefix = tempname;
	% 			copyfile(InFile, [NIIFileForARTPrefix '.nii.gz']);
	% 			gunzip([NIIFileForARTPrefix '.nii.gz']);
	% 			NIIFileForART = [NIIFileForARTPrefix '.nii'];
	% 		else
	% 			CompressedInput = false;
	% 			NIIFileForART = InFile;
	% 		end
	% 
	% 		NIIFileForART = fullfile(ARTTempDirectory, 'in_finalcrop.nii');
	% 		if(exist([NIIFileForART '.gz'], 'file') == 2)
	% 			gunzip([NIIFileForART '.gz']);
	% 			delete([NIIFileForART '.gz']);
	% 		end
	% 
	% 		NIIFileARTOutput = fullfile(ARTTempDirectory, 'out.nii');
	% 		ACPCDetect = '/usr/local/ART/bin/acpcdetect';
	% 		CommandString = ['export ARTHOME=/usr/local/ART; ' ACPCDetect ' -i ' NIIFileForART ' -o ' NIIFileARTOutput];
	% 		%disp(CommandString);
	% 		system(CommandString);

% 			NIIFileForART = fullfile(ARTTempDirectory, 'in_maxlabel_cropped.nii');
% 			if(exist([NIIFileForART '.gz'], 'file') == 2)
% 				gunzip([NIIFileForART '.gz']);
% 				delete([NIIFileForART '.gz']);
% 			end

			NIIFileARTOutput = fullfile(ARTTempDirectory, 'out_initial.nii');
			ACPCDetect = '/usr/local/ART/bin/acpcdetect';
			CommandString = ['export ARTHOME=/usr/local/ART; ' ACPCDetect ' -i ' NIIFileForART ' -o ' NIIFileARTOutput];
			%disp(CommandString);
			system(CommandString);
			%keyboard;
			%disp(fullfile(ARTTempDirectory, 'out.nii'))
	% 		if(CompressedInput)
	% 			delete([NIIFileForARTPrefix '.nii.gz']);
	% 			delete([NIIFileForARTPrefix '.nii']);
	% 		end
			[NII, AVW] = load_nii(NIIFileARTOutput);
			AVWMidSag = extract_mid_slice(AVW);
			NIIPixdims = NII.hdr.dime.pixdim(3:4);
			save(OutputMAT, 'AVWMidSag', 'NIIPixdims');
			%delete(fullfile(ARTTempDirectory, 'out.nii'));
			rmdir(ARTTempDirectory, 's');
		end
	end
elseif(isstruct(FileNameOrStruct))
	AVW = FileNameOrStruct.AVW;
	if(ismatrix(AVW))
		AVWMidSag = double(AVW);
	elseif(ndims(AVW) == 3)
		AVWMidSag = extract_mid_slice(AVW);
	end
	NIIPixdims = FileNameOrStruct.NIIPixdims;
end
%return;

%if(exist(InFile, 'file') ~= 2)
%	error(['Input file ' InFile ' is not a regular file or does not exist']);
%end
%keyboard;
[~, ~, ~] = mkdir(OutputDir);

OutputMAT = fullfile(OutputDir, [OutputPrefix '.mat']);
OutputPNG = fullfile(OutputDir, [OutputPrefix '.png']);

%[NII, AVW] = load_nii(InFile);

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
[~, FornixProbAVW] = load_nii('all_fornix_prob.nii.gz');
[~, ProbAVW] = load_nii('all_cc_prob.nii.gz');
[TemplateNII, TemplateAVW] = load_nii('all_msp_mean.nii.gz');
%[SmoothNII, SmoothAVW] = load_nii('all_msp_smooth_areas.nii.gz');

%DilationANGLES = load('all_cc_dilation_angle');
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

% if ndims(AVW) == 3
% 	NIIPixdims = NII.hdr.dime.pixdim(3:4);
% 
% 	if(mod(size(AVW, 2), 2) == 0)
% 		XI = size(AVW, 2) / 2;
% 		XFrac = 0.5;
% 		AVWMidSag = AVW(size(AVW, 1):-1:1, XI, size(AVW, 3):-1:1) .* (1 - XFrac) + ...
% 		AVW(size(AVW, 1):-1:1, XI + 1, size(AVW, 3):-1:1) .* XFrac;
% 		clear XI XFrac MidsagSliceIDX;
% 	else
% 		XI = ceil(size(AVW, 2) / 2);
% 		AVWMidSag = AVW(size(AVW, 1):-1:1, XI, size(AVW, 3):-1:1);
% 		clear XI;
% 	end
% 	AVWMidSag = squeeze(AVWMidSag)';
% 	AVWMidSag = double(AVWMidSag);
% elseif ndims(AVW) == 2
% 	AVWMidSag = double(AVW);
% 	T = NII.hdr.dime.pixdim(2:4);
% 	NIIPixdims = T(NII.hdr.dime.dim(2:4) > 1);
% 	clear T;
% end

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
%MatchedFilterSmoothAVW = imcrop(SmoothAVW, R);
%MatchedFilterAngles = imcrop(DilationANGLES.ANGLE, R);

LeftLimit = max(min(Templatexx), min(AVWMidSagxx));
UpperLimit = max(min(Templateyy), min(AVWMidSagyy));
RightLimit = min(max(Templatexx), max(AVWMidSagxx));
LowerLimit = min(max(Templateyy), max(AVWMidSagyy));

[ResampleX, ResampleY] = meshgrid(LeftLimit:TemplatePixdims(1):RightLimit, UpperLimit:TemplatePixdims(2):LowerLimit);

ResampledAVW = interp2(AVWMidSagX, AVWMidSagY, AVWMidSag, ResampleX, ResampleY);
if(GroundTruthGiven)
	ResampledGroundAVW = interp2(AVWMidSagX, AVWMidSagY, GroundTruthAVW, ResampleX, ResampleY, 'nearest');
else
	ResampledGroundAVW = [];
end

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
[cc2RRegMaxIDX] = find(cc2RRegMax);
Rofcc2RRegMax = R(cc2RRegMax);
[~, I] = min(Rofcc2RRegMax);

%cc2RRegMaxIMaxI = cc2RRegMaxI(I);
%cc2RRegMaxIMaxJ = cc2RRegMaxJ(I);

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
disp(['Expected location is of CC rank: ' num2str(ExpectedCCRank)]);
clear I;
%keyboard;

%FirstOffsetValid2R = find(yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1) & xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2), 1, 'first');
OffsetValid2R = find(...
	yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1) & ...
	xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2));

%total_offset2R = [yoffset2R(FirstOffsetValid2R), xoffset2R(FirstOffsetValid2R)];
%corr_offset2R = total_offset2R;
%xpeak2R = xoffset2R(FirstOffsetValid2R);
%ypeak2R = yoffset2R(FirstOffsetValid2R);

%real_total_offset = [yoffset2RReg, xoffset2RReg];
%real_total_offset = [ypeak2R, xpeak2R];
%realxpeak = xoffset2RReg;
%realypeak = yoffset2RReg;
%keyboard;

for MainFittingIteration = 1:2
	
	% if we are in the second fitting iteration then we are not using LK
	if MainFittingIteration == 1
		DoLK = true;
	else
		DoLK = false;
	end
	
	NumTries = 2;
	Offsets = zeros(NumTries + 1, 2);
	Offsets(1, :) = [yoffset2RReg, xoffset2RReg];
	Offsets(2:NumTries + 1, :) = [yoffset2R(OffsetValid2R(1:NumTries)), xoffset2R(OffsetValid2R(1:NumTries))];

	InitialResampledAVW = cell(1, NumTries + 1);
	MatchedFilterRemapped = cell(1, NumTries + 1);
	LKParameters = cell(1, NumTries + 1);
	LKCost = zeros(1, NumTries + 1);
	TX = cell(1, NumTries + 1);
	TY = cell(1, NumTries + 1);
	InterpX = cell(1, NumTries + 1);
	InterpY = cell(1, NumTries + 1);
	TemplateLKIMG = cell(1, NumTries + 1);
	TemplateProbLKIMG = cell(1, NumTries + 1);
	FornixProbLKIMG = cell(1, NumTries + 1);
	ResampledGroundCropped = cell(1, NumTries + 1);
	ResampledAVWCropped = cell(1, NumTries + 1);
	TemplateProbLKIMGCropped = cell(1, NumTries + 1);
	FornixProbLKIMGCropped = cell(1, NumTries + 1);
	OriginalOtsuMask = cell(1, NumTries + 1);
	OtsuMask = cell(1, NumTries + 1);
	TemplateOverlap = cell(1, NumTries + 1);

	for z = 1:NumTries + 1
		[InitialResampledAVW{z}, ...
		MatchedFilterRemapped{z}, ...
		LKParameters{z}, LKCost(z), ...
		TX{z}, TY{z}, InterpX{z}, InterpY{z}, ...
		TemplateLKIMG{z}, TemplateProbLKIMG{z}, FornixProbLKIMG{z}, ...
		ResampledGroundCropped{z}, ResampledAVWCropped{z}, ...
		TemplateProbLKIMGCropped{z}, FornixProbLKIMGCropped{z}, ...
		OriginalOtsuMask{z}, OtsuMask{z}, ...
		TemplateOverlap{z}] = do_lk_and_otsu(ResampledAVW, MatchedFilter, MatchedFilterProb, MatchedFilterFornixProb, ResampledGroundAVW, Offsets(z, :), DoLK);
	end
	% so sometimes the closest to the expected center fails and the centre needs to be the maximum correlation coefficient,
	% so new strategy, try the closest and if this fails (OtsuMask is empty) then resort to the CC values
	% find out which choice gives the best overlap on the template probabilities
	%%
	PCoverage = zeros(1, length(TemplateOverlap));
	%Props = {'Eccentricity', 'EquivDiameter', 'Area', 'Solidity'};
	DiceCoefficients = zeros(1, length(TemplateOverlap));
% 	HausdorffDistances = zeros(1, length(TemplateOverlap));
% 	DeformationNorm = zeros(1, length(TemplateOverlap));
% 	TransformationMatrix = cell(1, length(TemplateOverlap));
% 	ScalingCoefficients = zeros(2, length(TemplateOverlap));
% 	RotationMatrices = cell(1, length(TemplateOverlap));
% 	RotationAngles = zeros(1, length(TemplateOverlap));
	Distances = cell(1, length(TemplateOverlap));
	Angles = cell(1, length(TemplateOverlap));
	SumDistances = zeros(1, length(TemplateOverlap));
	SumAngles = zeros(1, length(TemplateOverlap));
	TransformAngleX = zeros(1, length(TemplateOverlap));
	TransformAngleY = zeros(1, length(TemplateOverlap));
	OtsuWithoutFornix = cell(1, length(TemplateOverlap));
	%keyboard;
	for z = 1:length(TemplateOverlap)
		PCoverage(z) = sum(TemplateProbLKIMGCropped{z}(:) .* OtsuMask{z}(:)) ./ sum(TemplateProbLKIMGCropped{z}(:));

% 		TransformationMatrix{z} = [1 + LKParameters{z}(1), LKParameters{z}(3), LKParameters{z}(5); ...
% 		LKParameters{z}(2), 1 + LKParameters{z}(4), LKParameters{z}(6); ...
% 		0, 0, 1];

% 		ScalingCoefficients(:, z) = sqrt(sum(TransformationMatrix{z}(1:2, 1:2) .* TransformationMatrix{z}(1:2, 1:2), 2));
% 		RotationMatrices{z} = bsxfun(@rdivide, TransformationMatrix{z}(1:2, 1:2), ScalingCoefficients(:, z));
% 		%keyboard;
% 		RotationAngles(z) = atan2(RotationMatrices{z}(2, 1), RotationMatrices{z}(2, 2));
		%TemplateWithoutFornix = TemplateProbLKIMGCropped{z} > 0.3;
		%TemplateProps = regionprops(, Props{:});
		%OtsuMaskProps = regionprops(double(OtsuMask{z}), Props{:});

		%for FIDX = 1:length(Props)
		%	MaskProps(z).(Props{FIDX}) = [TemplateProps.(Props{FIDX}) OtsuMaskProps.(Props{FIDX})];
		%end
		%clear TemplateProps OtsuMaskProps;
		T = TemplateProbLKIMGCropped{z} > 0.3;
		OtsuWithoutFornix{z} = OtsuMask{z} & ~FornixProbLKIMGCropped{z} > 0.3;
		OtsuWithoutFornix{z} = bwareaopen(OtsuWithoutFornix{z}, 50);
		DiceCoefficients(z) = dices_coefficient(T, OtsuWithoutFornix{z});
% 		HausdorffDistances(z) = hausdorff_distance(T, OtsuWithoutFornix{z});
% 		DeformationNorm(z) = norm(LKParameters{z}) + norm(LKParameters{z}([6, 5])' - Offsets(z, :));
		[Distances{z}, Angles{z}] = nearest_angle_distance(T, OtsuWithoutFornix{z});
		SumDistances(z) = sum(Distances{z});
		SumAngles(z) = sum(abs(Angles{z}));

		TransformAngleX(z) = atan2(TY{z}(1, 2) - TY{z}(1, 1), TX{z}(1, 2) - TX{z}(1, 1));
		TransformAngleY(z) = atan2(TY{z}(2, 1) - TY{z}(1, 1), TX{z}(2, 1) - TX{z}(1, 1)) - pi / 2;
		clear T;
	end

	%keyboard;
	%%
	% clf;
	% for z = 1:length(OtsuMask)
	% 	subplot(2, 3, z);
	% 	imshow(OtsuWithoutFornix{z}, []);
	% 	hold on;
	% 	[~, CC] = contour(TemplateProbLKIMGCropped{z}, [0.3, 0.3]);
	% 	set(CC, 'Color', 'r');
	% end
	% 
	%%
	% subplot 121;
	% imshow(OtsuMask{1}, []);
	% subplot 122;
	% imshow(OtsuMask{2}, []);
	%PCoverage
	%keyboard;
	ValidIDX = (TransformAngleX < pi / 3) & (TransformAngleY < pi / 3);
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
	[MinRigidFactor, MinRigidFactorI] = min(abs(TransformAngleX(ValidIDX) - TransformAngleY(ValidIDX)));
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

	% TestValues.Dice = MaxDiceCoefficientI;
	% Test
	% VS(1).
	% if(MaxDiceCoefficientI == MinHausdorffDistanceI)
	% 	BestLKI = MaxDiceCoefficientI;
	% else
	% 	keyboard;
	% end
	%keyboard;
	%BestLKI = ;
	DoCCGraphics = false;
	if(DoCCGraphics)
		clf;
		SR = 2;
		SC = NumTries + 1;
		%keyboard;
		%clf;
		for z = 1:NumTries + 1

		%	title({'Normxcorr result, red based on maximum cc,', 'blue based on minimum distance from expected'});
			[~, LOC] = ismember(z, ValidIDXIDX);

			subplot(SR, SC, z);
			imshow(ResampledAVW, []);
			hold on;
			SZ = size(MatchedFilter);

			T = line([Offsets(z, 2), Offsets(z, 2) + SZ(2)], [Offsets(z, 1), Offsets(z, 1)]); set(T, 'Color', 'b');
			T = line([Offsets(z, 2), Offsets(z, 2) + SZ(2)], [Offsets(z, 1) + SZ(1), Offsets(z, 1) + SZ(1)]); set(T, 'Color', 'b');
			T = line([Offsets(z, 2), Offsets(z, 2)], [Offsets(z, 1), Offsets(z, 1) + SZ(1)]); set(T, 'Color', 'b');
			T = line([Offsets(z, 2) + SZ(2), Offsets(z, 2) + SZ(2)], [Offsets(z, 1), Offsets(z, 1) + SZ(1)]); set(T, 'Color', 'b');

			T = line([TX{z}(1, 1) TX{z}(end, 1)], [TY{z}(1, 1), TY{z}(end, 1)]); set(T, 'Color', 'r');
			T = line([TX{z}(1, end) TX{z}(end, end)], [TY{z}(1, end), TY{z}(end, end)]); set(T, 'Color', 'r');
			T = line([TX{z}(1, 1) TX{z}(1, end)], [TY{z}(1, 1), TY{z}(1, end)]); set(T, 'Color', 'r');
			T = line([TX{z}(end, 1) TX{z}(end, end)], [TY{z}(end, 1), TY{z}(end, end)]); set(T, 'Color', 'r');
			title('Initial XCORR and LK alignment');
			%[MaxPCoverage, MaxPCoverageI] = max(PCoverage);
			%[MaxDiceCoefficient, MaxDiceCoefficientI] = max(DiceCoefficients);
			%[MinHausdorffDistance, MinHausdorffDistanceI] = min(HausdorffDistances);
			if(ValidIDX(z) == false)
				C = 'Invalid';
				TextColor = 'm';
			else
				C = cell(1, length(TestValues));
				for TestIDX = 1:length(TestValues)
					C{TestIDX} = [TestValues(TestIDX).text ': ' num2str(TestValues(TestIDX).values(z))];
					if(TestValues(TestIDX).bestidx == LOC)
						C{TestIDX} = ['{\bf' C{TestIDX} '}'];
					end
				end
				if BestAlignmentIDX == LOC
					TextColor = 'r';
				else
					TextColor = 'k';
				end
				if ~any(ValidIDX)
					TextColor = 'b';
				end
			end
			text(0, -0.01, C, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Color', TextColor);
			clear C;
			subplot(SR, SC, z + NumTries + 1);
			imshow(OriginalOtsuMask{z}, []);
			hold on;
			[~, CC] = contour(TemplateProbLKIMGCropped{z}, [0.3, 0.3]);
			%[~, CC] = contour(NewPSI, [0.3, 0.3]);
			set(CC, 'Color', 'r');
			xlabel('Otsu mask, template is red contour');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']

		end
	
		if MainFittingIteration == 1
			OutputFile = fullfile(OutputDir, 'cc', [OutputPrefix '_cc_choice.png']);
		else
			OutputFile = fullfile(OutputDir, 'cc', [OutputPrefix '_cc_fallback.png']);
		end
		[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
		%keyboard;
		FigPos = fullscreen_fig_pos;
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
		exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

		IMG = imread(OutputFile);
		IMG = imautocropwhite(IMG);
		imwrite(IMG, OutputFile);
	end
	%clf;
	%keyboard;
	% need to detect the situation where the LK transform has failed in all cases

	if(~all(isinf(SumDistances)) && any(ValidIDX))
		%keyboard;
		%return;
		%%%% THE HAUSDORFF DISTANCE DOESNT ALWAYS WORK %%%%
		%%%% Try and get the rotation angles of the LKParameters, choose the one that rotates the least
		%BestAlignmentIDX = ValidIDXIDX(BestAlignmentIDX)
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
		TemplateOverlap = TemplateOverlap{BestLKI};
		break;
	else
		if MainFittingIteration == 2
			disp('Something went wrong, check your file');
			return;
		end
	end
end

%keyboard;



%%
% clf;
% subplot 321;
% imshow(ResampledAVWCropped, []);
% subplot 323;
% imshow(ResampledAVWCropped, []);
% hold on;
% [I, J] = find(OtsuMask & (TemplateProbLKIMGCropped > 0.5));
% plot(J, I, '*');
% subplot(3, 2, 2);
% imshow(OtsuMask, []);
% subplot 325;
% imshow(TemplateProbLKIMGCropped, []);
% hold on;
% [I, J] = find(OtsuMask);
% plot(J, I, '*');

% Z = ResampledAVWCropped .* OtsuMask;
% Z(~OtsuMask) = NaN;
% 
% ZC = colfilt(Z, [25, 25], 'sliding', @mean_non_nan);
% clf;
% SR = 2;
% SC = 3;
% T = ResampledAVWCropped * 0;
% T(OtsuMask) = ResampledAVWCropped(OtsuMask) ./ ZC(OtsuMask);
% subplot(SR, SC, 1);
% imshow(T, []);
% subplot(SR, SC, 2);
% histfit(T(OtsuMask));
% subplot(SR, SC, 3);
% histfit(ResampledAVWCropped(OtsuMask));
% subplot(SR, SC, 4);
% imshow(T > 1.1, []);
% subplot(SR, SC, 5);
% imshow(ResampledAVWCropped, []);
% subplot(SR, SC, 6);
% [I, J] = find(T > 1.1);
% imshow(ResampledAVWCropped, []);
% hold on;
% plot(J, I, '*');
% %keyboard;
% 
% [~, ~, ~] = mkdir(fullfile(OutputDir, 'lowpass'));
% 
% OutputFile = fullfile(OutputDir, 'lowpass', [OutputPrefix '_lowpass_hist.png']);
% %keyboard;
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);
% return;

FirstStartMask = OtsuMask & (TemplateProbLKIMGCropped > 0.5);
%T = imerode(FirstStartMask, strel('disk', 1));
MeanIntensityInStartMask = mean(ResampledAVWCropped(FirstStartMask));
VarIntensityInStartMask = var(ResampledAVWCropped(FirstStartMask));

GaussianProb = normpdf(ResampledAVWCropped, MeanIntensityInStartMask, sqrt(VarIntensityInStartMask));
GaussianProb = GaussianProb ./ max(GaussianProb(:));
%GaussianProb(ResampledAVWCropped > MeanIntensityInStartMask) = 1;

LowPassMask = bwareaopen(OriginalOtsuMask, 200);
T = ResampledAVWCropped;
T(~LowPassMask) = NaN;
ZC = colfilt(T, [15, 15], 'sliding', @mean_non_nan);

ResampledAVWCroppedLowPass = zeros(size(ResampledAVWCropped));
ResampledAVWCroppedLowPass(LowPassMask) = ResampledAVWCropped(LowPassMask) ./ ZC(LowPassMask);
T = imerode(FirstStartMask, strel('disk', 1));
MeanIntensityInStartMaskLowPass = mean(ResampledAVWCroppedLowPass(T));
VarIntensityInStartMaskLowPass = var(ResampledAVWCroppedLowPass(T));

GaussianProbLowPass = normpdf(ResampledAVWCroppedLowPass, MeanIntensityInStartMaskLowPass, sqrt(VarIntensityInStartMaskLowPass));
GaussianProbLowPass = GaussianProbLowPass ./ max(GaussianProbLowPass(:));
GaussianProb(ResampledAVWCropped > MeanIntensityInStartMask) = 1;
clear T;

DoLowPassGraphics = false;
if(DoLowPassGraphics)
	
	clf;
	SR = 4;
	SC = 2;
	subplot(SR, SC, 1);
	M = [min(ResampledAVWCropped(FirstStartMask)) max(ResampledAVWCropped(FirstStartMask))];
	imshow(ResampledAVWCropped, M);
	subplot(SR, SC, 2);
	M = [min(ResampledAVWCroppedLowPass(FirstStartMask)) max(ResampledAVWCroppedLowPass(FirstStartMask))];
	imshow(ResampledAVWCroppedLowPass, M);
	subplot(SR, SC, 3);
	RGB = zeros([size(GaussianProbLowPass), 3]);
	T = GaussianProb .* (ResampledAVWCropped < MeanIntensityInStartMask);
	RGB(:, :, 1) = T;
	T = GaussianProb .* (ResampledAVWCropped >= MeanIntensityInStartMask);
	RGB(:, :, 3) = T;
	imshow(RGB, []);
	imshow(GaussianProb, []);
	subplot(SR, SC, 4);
	RGB = zeros([size(GaussianProbLowPass), 3]);
	T = GaussianProbLowPass .* (ResampledAVWCroppedLowPass < MeanIntensityInStartMaskLowPass);
	RGB(:, :, 1) = T;
	T = GaussianProbLowPass .* (ResampledAVWCroppedLowPass >= MeanIntensityInStartMaskLowPass);
	RGB(:, :, 3) = T;
	%imshow(RGB, []);
	imshow(GaussianProbLowPass, []);
	subplot(SR, SC, 5);
	histfit(ResampledAVWCropped(FirstStartMask));
	subplot(SR, SC, 6);
	histfit(ResampledAVWCroppedLowPass(FirstStartMask));
	OutputFile = fullfile(OutputDir, 'lowpass', [OutputPrefix '.png']);
	
 	[H, Bins] = hist(ResampledAVWCroppedLowPass(FirstStartMask), 100);
	F = gaussian_filter_max_1d(1);
	
	HD = imdilate(H(:), ones(7, 1));
	
 	Peaks = (HD == H(:));
% 	% find the highest peak
 	PeaksIDX = find(Peaks);
 	[PeakProb, I] = max(H(Peaks));
	MU = Bins(PeaksIDX(I));
	HF = conv(H(:), F(:), 'same');
	binwidth = Bins(2) - Bins(1);
	
	I = find(HF > PeakProb / 2, 1, 'last');
	IT = HF(I) - HF(I + 1);
	IFrac = (HF(I) - PeakProb / 2) / IT;
	FWHM = (Bins(I) + IFrac * binwidth - MU) * 2;
	SIGMA = FWHM ./ (2 * sqrt(2 * log(2)));
	%keyboard;
% 	init = [Bins(PeaksIDX(I)) / 2, Bins(PeaksIDX(I))];
% 	[label, model, llh] = emgm(ResampledAVWCroppedLowPass(FirstStartMask)', init);
 	subplot(SR, SC, 7);
 	bar(Bins, H);
 	hold on;
 	
 	area = sum(FirstStartMask(:)) * binwidth;
% 	y1 = normpdf(Bins, model.mu(1), model.Sigma(1)) * area;
% 	y2 = normpdf(Bins, model.mu(2), model.Sigma(2)) * area;
% 	plot(Bins, y1, Bins, y2);
% 	%keyboard;
	%y1 = normpdf(Bins, MeanIntensityInStartMaskLowPass, sqrt(VarIntensityInStartMaskLowPass)) * area;
	y1 = normpdf(Bins, MU, SIGMA) * area;
	plot(Bins, y1);
	subplot(SR, SC, 8);
	imshow(normpdf(ResampledAVWCroppedLowPass, MU, SIGMA), []);
	[~, ~, ~] = mkdir(fullfile(OutputDir, 'lowpass'));
	%keyboard;
	FigPos = fullscreen_fig_pos;
	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
	exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

	IMG = imread(OutputFile);
	IMG = imautocropwhite(IMG);
	imwrite(IMG, OutputFile);
end
%return;
% ButterWorthLowerCutoff = 3;
% ButterWorthHigherCutoff = ButterWorthLowerCutoff + min(size(ResampledAVWCropped)) / 4;
% ResampledAVWCroppedButter = butterworthbpf(ResampledAVWCropped, ButterWorthLowerCutoff, ButterWorthHigherCutoff, 2);
% 
% MeanIntensityInStartMaskButter = mean(ResampledAVWCroppedButter(FirstStartMask));
% VarIntensityInStartMaskButter = var(ResampledAVWCroppedButter(FirstStartMask));
% 
% GaussianProbButter = normpdf(ResampledAVWCroppedButter, MeanIntensityInStartMaskButter, sqrt(VarIntensityInStartMaskButter));
% GaussianProbButter = GaussianProbButter ./ max(GaussianProbButter(:));
% %GaussianProb(ResampledAVWCropped > MeanIntensityInStartMask) = 1;

% clf;
% SR = 2;
% SC = 3;
% subplot(SR, SC, 1);
% imshow(ResampledAVWCropped, []);
% hold on;
% [C, CC] = contour(FirstStartMask, [0.5, 0.5]);
% set(CC, 'Color', 'r');
% subplot(SR, SC, 2);
% imshow(ResampledAVWCroppedButter, []);
% subplot(SR, SC, 3);
% imshow(GaussianProb, []);
% subplot(SR, SC, 4);
% imshow(GaussianProbButter, []);
% subplot(SR, SC, 5);
% histfit(ResampledAVWCropped(FirstStartMask));
% subplot(SR, SC, 6);
% histfit(ResampledAVWCroppedButter(FirstStartMask));
% 
% 
% [~, ~, ~] = mkdir(fullfile(OutputDir, 'butter'));
% 
% OutputFile = fullfile(OutputDir, 'butter', [OutputPrefix '_prob_butter.png']);
% %keyboard;
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);
% 
% return;

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
%ExclusionWeight = 
ExclusionA = min(FornixProbLKIMGCropped * 50, 1) .* GaussianProb;
ExclusionB = GaussianGradMAG;
ExclusionC = exp(-2 .* GaussianProb .* (TemplateProbLKIMGCropped + 2));
ExclusionA = ExclusionA ./ max(ExclusionA(:));
ExclusionB = ExclusionB ./ max(ExclusionB(:));
ExclusionC = ExclusionC ./ max(ExclusionC(:));

%ExclusionWeight = min(FornixProbLKIMGCropped * 50, 1) .* GaussianProb + GaussianGradMAG + exp(-2*GaussianProb);
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
clear StrongI StrongJ;

%L = regionprops(StartMask, 'Area');
%StartMaskAreas = [L.Area];

StartMask = bwareaopen(StartMask, 200);
%subplot 122;
% imshow(StartMask, []);
OriginalPSI = double(StartMask);
%OriginalS = OtsuMask & (TemplateProbLKIMGCropped > 0.2) & (GaussianProb > 0.3 | ResampledAVWCropped > MeanIntensityInStartMask);
%StartMask;

%T = S + FornixProbLKIMGCropped;
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

	T = bw_join_mst_by_boundaries(StartMask);
	StartMask = StartMask | (T & GaussianProb > 0.3);
	DMask = DMask | (imdilate(T, strel('disk', 2)) & GaussianProb > 0.2);
	clear T;

end

DMaskCC = bwconncomp(DMask);
DMaskL = labelmatrix(DMaskCC);
% 	imshow(label2rgb(DMaskL, 'lines', 'k'), []);
% 	
% find the region in DMask that mostly overlaps with the startmask, filter out other regions
DMaskR = regionprops(DMaskCC, 'PixelIdxList');

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
%	keyboard;
	%imshow(label2rgb(DMaskL, 'lines', 'k'), []);
	disp('overlaps dont match, the DMask needs to be joined');
	NewDMask = bwareaopen(DMask, 200) & ismember(DMaskL, find(Overlaps > 0.1));
	return;
	NewDMaskJoiningSegments = bw_join_mst_by_boundaries(NewDMask);
	NewDMaskJoined = NewDMask | NewDMaskJoiningSegments;
	DMask = NewDMaskJoined;
	%keyboard;
else
	[~, L] = max(Overlaps);
	DMask = (DMaskL == L);
end
clear Overlaps DMaskL L;

SE = strel('disk', 4);
DMask = imclose(DMask, SE);
StartMask = imclose(StartMask, SE);

MeanIntensityInFinalStartMask = mean(ResampledAVWCropped(StartMask));
VarIntensityInFinalStartMask = var(ResampledAVWCropped(StartMask));

GaussianProbFinal = normpdf(ResampledAVWCropped, MeanIntensityInFinalStartMask, sqrt(VarIntensityInFinalStartMask));
GaussianProbFinal = GaussianProbFinal ./ max(GaussianProbFinal(:));
%GaussianProbFinalNoPositive = GaussianProbFinal;
GaussianProbFinal(ResampledAVWCropped > MeanIntensityInFinalStartMask) = 1;


StartMask = StartMask & (GaussianProbFinal > 0.01);
DMask = DMask & (GaussianProbFinal > 0.001);

DMaskCC = bwconncomp(DMask);
DMaskL = labelmatrix(DMaskCC);
DMaskAreas = cellfun('length', DMaskCC.PixelIdxList);
[~, L] = max(DMaskAreas);
%keyboard;
DMask = (DMaskL == L);

StartMaskCC = bwconncomp(StartMask);
StartMaskL = labelmatrix(StartMaskCC);
StartMaskAreas = cellfun('length', StartMaskCC.PixelIdxList);
[~, L] = max(StartMaskAreas);
StartMask = (StartMaskL == L);

DMask = bwfill(DMask, 'holes');
StartMask = bwfill(StartMask, 'holes');

% [~, Angles, SmoothedBoundaryA, SmoothedBoundaryB] = nearest_angle_distance(DMask, TemplateProbLKIMGCropped > 0.3);
% 
% OutlierPixels = false(size(DMask));
% 
% for z = 1:length(SmoothedBoundaryA)
% 	if(abs(Angles(z)) > pi / 4)
% 		OutlierPixels(floor(SmoothedBoundaryA(z, 1)), floor(SmoothedBoundaryA(z, 2))) = 1;
% 		OutlierPixels(floor(SmoothedBoundaryA(z, 1)), ceil(SmoothedBoundaryA(z, 2))) = 1;
% 		OutlierPixels(ceil(SmoothedBoundaryA(z, 1)), floor(SmoothedBoundaryA(z, 2))) = 1;
% 		OutlierPixels(ceil(SmoothedBoundaryA(z, 1)), ceil(SmoothedBoundaryA(z, 2))) = 1;
% 	end
% end
% %keyboard;
% clf;
% SR = 2;
% SC = 2;
% subplot(SR, SC, 1);
% imshow(ResampledAVWCropped, []);
% hold on;
% [I, J] = find(OutlierPixels);
% plot(J, I, '*');
% S = SmoothedBoundaryA(:, [2 1]);
% SS = streamline({S});
% set(SS, 'Color', 'r');
% subplot(SR, SC, 2);
% imshow(DMask, []);
% subplot(SR, SC, 3);
% histfit(ResampledAVWCropped(DMask));
% subplot(SR, SC, 4);
% imshow(ResampledAVWCropped, [min(ResampledAVWCropped(DMask)), max(ResampledAVWCropped(DMask))]);
% [~, ~, ~] = mkdir(fullfile(OutputDir, 'outliers'));
% 
% OutputFile = fullfile(OutputDir, 'outliers', [OutputPrefix '_dmask_outliers.png']);
% %keyboard;
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);
% 
% return;
%keyboard;
% 
% subplot 221;
% imshow(GaussianProb, []);
% title(num2str([MeanIntensityInStartMask, VarIntensityInStartMask]));
% 
% subplot 222;
% imshow(GaussianProbFinal, []);
% title(num2str([MeanIntensityInFinalStartMask, VarIntensityInFinalStartMask]));
% hold on;
% [~, CC] = contour(DMask, [0.5 0.5]);
% set(CC, 'Color', 'r');
% [~, CC] = contour(StartMask, [0.5 0.5]);
% set(CC, 'Color', 'b');
% %[~, CC] = contour(GaussianProbFinal > 0.001, [0.5 0.5]);
% %set(CC, 'Color', 'b');
% 
% %%
% subplot 223;
% imshow(GaussianProbFinal > 0.001, []);
%%

%keyboard;
OriginalPSI = zeros(size(StartMask));
OriginalPSI(StartMask) = 1;
OriginalPSI(~StartMask) = -1;
%OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);
%OriginalPSI = computeDistanceFunction2d(OriginalPSI, [1, 1], [], 2);

% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% return;

% D = double(DMask);
% D(DMask) = 1;
% D(~DMask) = -1;
%D = computeDistanceFunction2d(D, [1, 1], [], 2);
%clear S;

%keyboard;
%FornixProbLKIMGCropped(isinf(FornixProbLKIMGCropped)) = 0;
%W = 5;
%[NewPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(-D, OriginalPSI, W - (TemplateProbLKIMGCropped * (W - 1)), double(ExclusionWeight), [1, 1], 1000, 200, 10, 5, ResampledAVWCropped);
%[NewPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(-D, OriginalPSI, ones(size(D)), double(ExclusionWeight), [1, 1], 1000, 200, 10, 5, ResampledAVWCropped);
%[NewPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(-D, OriginalPSI, ones(size(D)), zeros(size(D)), [1, 1], 1000, 200, 10, 5, ResampledAVWCropped);
%[NewPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(-D, OriginalPSI, ones(size(D)), zeros(size(D)), [1, 1], 1000, 200, 5, 5, ResampledAVWCropped);
% no illusory contours
NewPSI = OriginalPSI;
FinalSeg = StartMask;%NewPSI > 0;
% T = ResampledAVWCropped(FinalSeg);
% MeanIntensityInSEG = mean(T);
% VarIntensityInSEG = var(T);
% clear T;
% % F = F(1:round(0.95 * length(F)));
% % 
% % MeanIntensityInSEGR = mean(F);
% % VarIntensityInSEGR = var(F);
% % 
% % GaussianProbSEGR = normpdf(ResampledAVWCropped, MeanIntensityInSEGR, sqrt(VarIntensityInSEGR));
% % %GaussianProbSEG = GaussianProbSEG ./ max(GaussianProbSEG(:));
% % %GaussianProbSEG(ResampledAVWCropped > MeanIntensityInSEG) = 1;
% % GaussianProbSEGR(~FinalSeg) = 0;
% 
% GaussianProbSEG = normpdf(ResampledAVWCropped, MeanIntensityInSEG, sqrt(VarIntensityInSEG));
% GaussianProbSEG = GaussianProbSEG ./ max(GaussianProbSEG(:));
% GaussianProbSEG(ResampledAVWCropped < MeanIntensityInSEG) = 1;
% GaussianProbSEG(~FinalSeg) = 1;
% %%
% % subplot 221;
% % imshow(ResampledAVWCropped, []);
% % subplot 222;
% % imshow(GaussianProbFinal, []);
% % subplot 223;
% % imshow(GaussianProbSEG > 0.25, []);
% 
% FinalSegOpened = imopen(GaussianProbSEG > 0.05, strel('disk', 2));
% FinalSegOpened(~FinalSeg) = 0;
% FinalSegOpenedCC = bwconncomp(FinalSegOpened);
% FinalSegOpenedAreas = cellfun('length', FinalSegOpenedCC.PixelIdxList);
% FinalSegOpenedL = labelmatrix(FinalSegOpenedCC);
% [~, I] = max(FinalSegOpenedAreas);
% FinalSegOpenedM = FinalSegOpenedL == I;
% FinalSegOpenedM = bwfill(FinalSegOpenedM, 'holes');
% subplot 224;
% imshow(ResampledAVWCropped, []);
% hold on;
% [~, CC] = contour(FinalSegOpenedM, [0.5, 0.5]);
% set(CC, 'Color', 'r');
%imshow(medfilt2(GaussianProbSEGR, [3, 3]), []);
%hist(ResampledAVWCropped(FinalSeg), 500);
%%

%%
% tophat - bottomhat test

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
%%
% subplot 221;
% imshow(TopHat, []);
% subplot 222;
% imshow(BottomHat, []);

%T = BottomHat + TopHat;
T = BottomHat + TopHat;
BottomTopHat = BottomHat + TopHat;
%T(:) = 0;
% subplot 223;
% imshow(T, []);
% subplot 224;
%M = median(T(FinalSeg));
% imshow(T > M * 3, []);

%FinalSegFirstCull = FinalSeg & (T < M * 3);
%clear T;
%S = ResampledAVWCropped > MeanIntensityInFinalStartMask + 0.5 * sqrt(VarIntensityInFinalStartMask);
%%
OutlierPixelsStrong = T > (median(T(FinalSeg)) * 3);
OutlierPixelsWeak = T > (median(T(FinalSeg)) * 2);
% 
% [OutlierPixelsStrongR, OutlierPixelsStrongC] = find(OutlierPixelsStrong);
% OutlierPixels = bwselect(OutlierPixelsWeak, OutlierPixelsStrongC, OutlierPixelsStrongR, 8);
% clf;
% subplot 131;
% imshow(OutlierPixelsStrong, []);
% subplot 132;
% imshow(OutlierPixelsWeak, []);
% subplot 133;
% imshow(OutlierPixels, []);
% %%
% keyboard;

%%
% [~, Angles, SmoothedBoundaryA, SmoothedBoundaryB] = nearest_angle_distance(FinalSeg, TemplateProbLKIMGCropped > 0.3);
% hold on;
% CMAPX = linspace(0, pi, 256);
% CMAP = jet(256);
% 
% AngleRGB = zeros(numel(Angles), 3);
% AngleRGB(:, 1) = interp1(CMAPX, CMAP(:, 1), Angles);
% AngleRGB(:, 2) = interp1(CMAPX, CMAP(:, 2), Angles);
% AngleRGB(:, 3) = interp1(CMAPX, CMAP(:, 3), Angles);
% 
% for z = 1:length(SmoothedBoundaryA)
% 	IDX = z + 1;
% 	if(IDX > size(SmoothedBoundaryA, 1))
% 		IDX = 1;
% 	end
% 	L = line([SmoothedBoundaryA(z, 2) SmoothedBoundaryA(IDX, 2)], [SmoothedBoundaryA(z, 1) SmoothedBoundaryA(IDX, 1)]);
% 	set(L, 'Color', AngleRGB(z, :), 'LineWidth', 3);
% end
%OutlierPixels = false(size(OutlierPixelsWeak));
%SmoothedBoundaryB = SmoothedBoundaryB(:, [2 1]);
%S = streamline({SmoothedBoundaryB});
%set(S, 'Color', 'r');

% for z = 1:length(SmoothedBoundaryA)
% 	if(abs(Angles(z)) > pi / 4)
% 		OutlierPixels(floor(SmoothedBoundaryA(z, 1)), floor(SmoothedBoundaryA(z, 2))) = 1;
% 		OutlierPixels(floor(SmoothedBoundaryA(z, 1)), ceil(SmoothedBoundaryA(z, 2))) = 1;
% 		OutlierPixels(ceil(SmoothedBoundaryA(z, 1)), floor(SmoothedBoundaryA(z, 2))) = 1;
% 		OutlierPixels(ceil(SmoothedBoundaryA(z, 1)), ceil(SmoothedBoundaryA(z, 2))) = 1;
% 	end
% end

if(any(OutlierPixelsStrong(:)))
	[OutlierPixelsStrongR, OutlierPixelsStrongC] = find(OutlierPixelsStrong);
	OutlierPixels = bwselect(OutlierPixelsWeak, OutlierPixelsStrongC, OutlierPixelsStrongR, 8);
	%OutlierPixels = OutlierPixelsStrong;
% 	clf;
% 	subplot 221;
% 	imshow(OutlierPixelsStrong, []);
% 	subplot 222;
% 	imshow(OutlierPixelsWeak, []);
% 	subplot 223;
% 	imshow(OutlierPixels, []);
% 	subplot 224;
% 	imshow(OutlierPixelsStrongOld, []);
% 	keyboard;
else
	OutlierPixels = false(size(T));
end

% dilate the outlier pixels
if(~any(OutlierPixels(:)))
	FinalSegOpenedM	= FinalSeg;
else
	SE = strel('disk', 3);
	OutlierPixelsDilated = imdilate(OutlierPixels, SE);
	FinalSegWithout = FinalSeg & ~OutlierPixelsDilated;
	T = imdilate(OutlierPixelsDilated, ones(3, 3));
	FinalSegWithoutBorders = T & FinalSegWithout;
	
	FinalSegOpened = imdilate(FinalSegWithoutBorders, SE) | FinalSegWithout;%imopen(FinalSegFirstCull, strel('disk', 3));

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
	
% 	clf;
% 	subplot 121;
% 	imshow(FinalSegOpened, []);
% 	subplot 122;
% 	imshow(FinalSegOpenedG, []);
% 	keyboard;
	[~, I] = max(FinalSegOpenedAreas);
	FinalSegOpenedM = FinalSegOpenedL == I;
	FinalSegOpenedM = imclose(FinalSegOpenedM, SE);
	FinalSegOpenedM = bwfill(FinalSegOpenedM, 'holes');
	
end

FinalSegOpenedMOpened = imopen(FinalSegOpenedM, strel('disk', 1));
FinalSegOpenedMOpenedLost = ~FinalSegOpenedMOpened & FinalSegOpenedM;

FinalSegOpenedMOpenedLostCC = bwconncomp(FinalSegOpenedMOpenedLost, 8);
%FinalSegOpenedMOpenedLostL = labelmatrix(FinalSegOpenedMOpenedLostCC);
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
% %keyboard;
% 
% %%
% 
% %%
%clear FinalSegOpenedMOpened FinalSegOpenedMOpenedLostCC FinalSegOpenedMOpenedLost;
%OldFinalSegOpenedMThickened = FinalSegOpenedMThickened;
%OutlierPixels = bwareaopen(OutlierPixels, 6);
%OutlierPixels(:) = 0;%bwareaopen(OutlierPixels, 6);
FinalSegOpenedMThickened = imdilate(FinalSegOpenedM, strel('disk', 1));
FinalSegOpenedMThickened = bwfill(FinalSegOpenedMThickened, 'holes', 8);

FinalSeg = imdilate(FinalSeg, strel('disk', 1));
FinalSeg = bwfill(FinalSeg, 'holes', 8);

%%
% clf;
% SR = 4;
% SC = 2;
% %subplot(SR, SC, 1);
% imshow(ResampledAVWCropped, []);
% hold on;
% [I, J] = find(OutlierPixels);
% plot(J, I, '*');
% S = SmoothedBoundaryA(:, [2 1]);
% SS = streamline({S});
% set(SS, 'Color', 'r');
% [~, CC] = contour(FinalSegOpenedM, [0.5, 0.5]);
% set(CC, 'Color', 'r');
% [~, CC] = contour(FinalSegOpenedMOpened, [0.5, 0.5]);
% set(CC, 'Color', 'b');
% [~, CC] = contour(FinalSeg, [0.5, 0.5]);
% set(CC, 'Color', 'm');
% [~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
%set(CC, 'Color', 'g');
%subplot(SR, SC, 2);
%imshow(OutlierPixels, []);

% subplot(SR, SC, 3);
% imshow(BottomTopHat, []);
% subplot(SR, SC, 4);
% hist(BottomTopHat(BottomTopHat ~= 0), 100);
% subplot(SR, SC, 5);
% imshow(TopHat, []);
% subplot(SR, SC, 6);
% imshow(BottomHat, []);
% 
% subplot(SR, SC, 7);
% imshow(GaussianGradMAG > 0.2, []);
% %subplot(SR, SC, 8);
% %clf;
% imshow(GaussianProbFinalNoPositive > 0.2 & FinalSeg, []);

%%
% imshow(ResampledAVWCropped > MeanIntensityInFinalStartMask + 2 * sqrt(VarIntensityInFinalStartMask), []);
% %%


%%
%keyboard;
%%
% subplot 121;
% imshow(FinalSegOpenedM, []);
% subplot 122;
% imshow(FinalSegOpenedMThickened, []);
% keyboard;
%hist(TopHat(FinalSeg), 100);
%subplot 224;
%hist(BottomHat(FinalSeg), 100);



%G = wblfit(TopHat(FinalSeg));
%Y = wblpdf(TopHat(FinalSeg), G(1), G(2));
%T = double(FinalSeg);
%T(FinalSeg) = Y;
%imshow(T, []);
%keyboard;
DoThicknesses = false;

if(DoThicknesses)
	ThicknessCalculationError = false;
	try
		[Contours, SolvedImage, InitialPoints, ArcLengthsLeft, ArcLengthsRight, LeftOffsetIDXStart, CurLeftOffsetIDX, RightOffsetIDXStart, CurRightOffsetIDX] = endpoint_find_one_subject(FinalSegOpenedMThickened);

		[MaskInnerBoundary, ...
			MaskOuterBoundary, ...
			MaskFree, ...
			SolvedImage, ...
			XY, ...
			X, ...
			Y, ...
			NormFX, ...
			NormFY, ...
			ValidStreamlines, ...
			StartV] = laplace_get_points_2d_auto_mw(Contours.xi, Contours.yi, Contours.xo(end:-1:1), Contours.yo(end:-1:1), 1, 1, 1, 100);

		XYForArcLengths = cell(size(XY));

		for z = 1:length(XY)
			XYForArcLengths{z} = bsxfun(@times, XY{z}, TemplatePixdims(:)');
		end

		T = cellfun(@arc_length, XYForArcLengths(ValidStreamlines));
		ArcLengths = zeros(size(XY));
		ArcLengths(ValidStreamlines) = T;
		clear T;
	catch ME
		disp('There was a problem with the thickness calculations');
		%ME.stack
		ThicknessCalculationError = true;
		disp(ME.message);
		disp('Stack:');
		for z = 1:length(ME.stack)
			disp([ME.stack(z).file '>' ME.stack(z).name ': line ' num2str(ME.stack(z).line)]);
		end
	end

	if(GroundTruthGiven)
		%keyboard;
		try
			ThicknessCalculationErrorGround = false;
			[GND.Contours, GND.SolvedImage, GND.InitialPoints, GND.ArcLengthsLeft, GND.ArcLengthsRight, GND.LeftOffsetIDXStart, GND.CurLeftOffsetIDX, GND.RightOffsetIDXStart, GND.CurRightOffsetIDX] = endpoint_find_one_subject(imdilate(ResampledGroundCropped, strel('disk', 1)));
			%keyboard;
			[GND.MaskInnerBoundary, ...
				GND.MaskOuterBoundary, ...
				GND.MaskFree, ...
				GND.SolvedImage, ...
				GND.XY, ...
				GND.X, ...
				GND.Y, ...
				GND.NormFX, ...
				GND.NormFY, ...
				GND.ValidStreamlines, ...
				GND.StartV] = laplace_get_points_2d_auto_mw(GND.Contours.xi, GND.Contours.yi, GND.Contours.xo(end:-1:1), GND.Contours.yo(end:-1:1), 1, 1, 1, 100);
				XYForArcLengths = cell(size(XY));
			%keyboard;
			for z = 1:length(XY)
				XYForArcLengths{z} = bsxfun(@times, GND.XY{z}, TemplatePixdims(:)');
			end

			T = cellfun(@arc_length, XYForArcLengths(GND.ValidStreamlines));
			GND.ArcLengths = zeros(size(GND.XY));
			GND.ArcLengths(GND.ValidStreamlines) = T;
			clear T;
		catch ME
			disp('There was a problem with the thickness calculations, ground truth image');
			ThicknessCalculationErrorGround = true;
			disp(ME.message);
			disp('Stack:');
			for z = 1:length(ME.stack)
				disp([ME.stack(z).file '>' ME.stack(z).name ': line ' num2str(ME.stack(z).line)]);
			end
		end
	end
end

PlotGraphics = true;

if(PlotGraphics)
	
	%%
	clf;
	SR = 3;
	SC = 4;
	
	SZ = size(MatchedFilter);
	%keyboard;
% 	subplot(SR, SC, 1);
% 	imshow(cc2R, []);
% 	hold on;
	%plot(Offsets, SortedCCID(FirstOffsetValid2R), 'r*');
	%plot(cc2RRegMaxIMaxJ, cc2RRegMaxIMaxI, 'b*');
	
	
	%title({'Normxcorr result'});

	subplot(SR, SC, 1);
	imshow(cc2R, []);
	%[I, J] = find
	hold on;
	plot(cc2RRegMaxJ, cc2RRegMaxI, '*');
	plot(Offsets(BestLKI, 2) + SZ(2), Offsets(BestLKI, 1) + SZ(1), 'r*');
	title('Regional maxima of cc2R');
	
	subplot(SR, SC, 2);
	imshow(ResampledAVW, []);
	hold on;
	
	T = line([Offsets(BestLKI, 2), Offsets(BestLKI, 2) + SZ(2)], [Offsets(BestLKI, 1), Offsets(BestLKI, 1)]); set(T, 'Color', 'b');
	T = line([Offsets(BestLKI, 2), Offsets(BestLKI, 2) + SZ(2)], [Offsets(BestLKI, 1) + SZ(1), Offsets(BestLKI, 1) + SZ(1)]); set(T, 'Color', 'b');
	T = line([Offsets(BestLKI, 2), Offsets(BestLKI, 2)], [Offsets(BestLKI, 1), Offsets(BestLKI, 1) + SZ(1)]); set(T, 'Color', 'b');
	T = line([Offsets(BestLKI, 2) + SZ(2), Offsets(BestLKI, 2) + SZ(2)], [Offsets(BestLKI, 1), Offsets(BestLKI, 1) + SZ(1)]); set(T, 'Color', 'b');
	% T = line([Offsets(BestLKI, 2)2, Offsets(BestLKI, 2)2 + SZ(2)], [Offsets(BestLKI, 1)2, Offsets(BestLKI, 1)2]); set(T, 'Color', 'g');
	% T = line([Offsets(BestLKI, 2)2, Offsets(BestLKI, 2)2 + SZ(2)], [Offsets(BestLKI, 1)2 + SZ(1), Offsets(BestLKI, 1)2 + SZ(1)]); set(T, 'Color', 'g');
	% T = line([Offsets(BestLKI, 2)2, Offsets(BestLKI, 2)2], [Offsets(BestLKI, 1)2, Offsets(BestLKI, 1)2 + SZ(1)]); set(T, 'Color', 'g');
	% T = line([Offsets(BestLKI, 2)2 + SZ(2), Offsets(BestLKI, 2)2 + SZ(2)], [Offsets(BestLKI, 1)2, Offsets(BestLKI, 1)2 + SZ(1)]); set(T, 'Color', 'g');

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

	subplot(SR, SC, 3);
	imshow(OriginalOtsuMask, []);
	hold on;
	[~, CC] = contour(OriginalPSI, [0, 0]);
	set(CC, 'Color', 'b');
	title('Original Otsu');
	
	subplot(SR, SC, 4);
	imshow(ResampledAVWCropped, []);
	hold on;
	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
	set(CC, 'Color', 'r');
	%[~, CC] = contour(OriginalPSI, [0, 0]);
	%set(CC, 'Color', 'b');
	[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	set(CC, 'Color', 'g');
	title('Original Image');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']

	subplot(SR, SC, 5);
	imshow(GaussianProb, []);
	hold on;
	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
	set(CC, 'Color', 'r');
	%[~, CC] = contour(OriginalPSI, [0, 0]);
	%set(CC, 'Color', 'b');
	[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	set(CC, 'Color', 'g');
	title('Gaussian probability of pixels');

	subplot(SR, SC, 6);
	imshow(TemplateProbLKIMGCropped, []);
	hold on;
	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
	[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	set(CC, 'Color', 'g');
	title('Probability of CC from template');
	text(0.5, -0.01, {'In all these images', 'red: result', 'green: artefacts removed'}, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
% 	subplot(SR, SC, 7);
% 	imshow(TemplateProbLKIMGCropped, []);
% 	hold on;
% 	[~, CC] = contour(OriginalOtsuMask, [0.5, 0.5]);
% 	set(CC, 'Color', 'r');
% 	%title(['Probability of CC from template with original otsu borders, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']);
% 	title(['Probability of CC from template with original otsu borders']);

% 	subplot(SR, SC, 8);
% 	imshow(OriginalOtsuSeg, []);
% 	hold on;
% 	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
% 	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
% 	title(['OtsuMask ' num2str(OriginalOtsuSegCounts(1:3)', '%d ')]);

	subplot(SR, SC, 7);
	imshow(ExclusionWeight, []);
	hold on;
	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
	set(CC, 'Color', 'r');
	%[~, CC] = contour(OriginalPSI, [0, 0]);
	%set(CC, 'Color', 'b');
	[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	set(CC, 'Color', 'g');
	title('Exclusion mask');
%%	
	if(DoThicknesses)
		if(~ThicknessCalculationError)
			subplot(SR, SC, 8);
			plot(ArcLengths);
			xlabel('Node');
			ylabel('Thickness (mm)');
			title('Thickness');

			subplot(SR, SC, 9);
			imshow(ResampledAVWCropped, []);
			hold on;
			%keyboard;
			S = streamline(XY(ValidStreamlines));
			set(S, 'Color', 'r');
			plot(Contours.xo, Contours.yo, Contours.xi, Contours.yi);
			title('Streamlines');
		end
	end
	if(GroundTruthGiven)
		subplot(SR, SC, 10);
		imshow(ResampledAVWCropped, []);	
		hold on;
		[~, CC] = contour(ResampledGroundCropped, [0.5, 0.5]);
		set(CC, 'Color', 'b');
		[~, CC] = contour(FinalSeg, [0.5, 0.5]);
		set(CC, 'Color', 'r');
		[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
		set(CC, 'Color', 'g');
		title('Ground truth (blue) + Estimated (red) + Final (green)');		
		if(DoThicknesses)
			if(~ThicknessCalculationErrorGround)
				subplot(SR, SC, 11);
				imshow(ResampledAVWCropped, []);
				hold on;
				S = streamline(XY(ValidStreamlines));
				set(S, 'Color', 'r');
		%		keyboard;
				S = streamline(GND.XY(GND.ValidStreamlines));
				set(S, 'Color', 'b');
				title('Streamlines + Ground Truth');

				subplot(SR, SC, 12);
				plot([ArcLengths, GND.ArcLengths(:)]);
				xlabel('Node');
				ylabel('Thickness (mm)');
				title('Thickness + Ground Truth');
			end
		end
	end
% 	clf;
% 
% 	imshow(ResampledAVWCropped, []);
% 	hold on;
% 	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
% 	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
% 	[~, CC] = contour(FinalSegOpenedM, [0.5, 0.5]);
% 	set(CC, 'Color', 'g');
% 	title('Original Image');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
%	test for opening method	
% 	subplot 131;
% 	imshow(OutlierPixels, []);
% 	subplot 132;
% 	imshow(ResampledAVWCropped, []);
% 	hold on;
% 	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
% 	set(CC, 'Color', 'r');
% 	[~, CC] = contour(OriginalPSI, [0, 0]);
% 	set(CC, 'Color', 'b');
% 	[~, CC] = contour(FinalSegOpenedM, [0.5, 0.5]);
% 	set(CC, 'Color', 'g');
% 	title('Original Image');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
% 	subplot 133;
% 	imshow(FinalSegOpenedM, []);
	%keyboard;
	FigPos = fullscreen_fig_pos;
	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
	exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
	
	IMG = imread(OutputPNG);
	IMG = imautocropwhite(IMG);
	imwrite(IMG, OutputPNG);
	%delete(gcf);
end

OutputS.InitialSeg = OriginalPSI > 0;
OutputS.FinalSeg = FinalSeg;
OutputS.FinalSegArtefactsRemoved = FinalSegOpenedMThickened;
OutputS.IMG = ResampledAVWCropped;
%keyboard;
if(DoThicknesses)
	if(~ThicknessCalculationError)
		[OutputS.FinalThickness.Contours, OutputS.FinalThickness.SolvedImage, OutputS.FinalThickness.InitialPoints, OutputS.FinalThickness.ArcLengthsLeft, OutputS.FinalThickness.ArcLengthsRight, OutputS.FinalThickness.LeftOffsetIDXStart, OutputS.FinalThickness.CurLeftOffsetIDX, OutputS.FinalThickness.RightOffsetIDXStart, OutputS.FinalThickness.CurRightOffsetIDX] = ...
		deal(Contours, SolvedImage, InitialPoints, ArcLengthsLeft, ArcLengthsRight, LeftOffsetIDXStart, CurLeftOffsetIDX, RightOffsetIDXStart, CurRightOffsetIDX);

		[OutputS.FinalThickness.MaskInnerBoundary, ...
		OutputS.FinalThickness.MaskOuterBoundary, ...
		OutputS.FinalThickness.MaskFree, ...
		OutputS.FinalThickness.SolvedImage, ...
		OutputS.FinalThickness.XY, ...
		OutputS.FinalThickness.X, ...
		OutputS.FinalThickness.Y, ...
		OutputS.FinalThickness.NormFX, ...
		OutputS.FinalThickness.NormFY, ...
		OutputS.FinalThickness.ValidStreamlines, ...
		OutputS.FinalThickness.StartV] = deal(MaskInnerBoundary, ...
		MaskOuterBoundary, ...
		MaskFree, ...
		SolvedImage, ...
		XY, ...
		X, ...
		Y, ...
		NormFX, ...
		NormFY, ...
		ValidStreamlines, ...
		StartV);
		OutputS.FinalThicknesses = ArcLengths;
	end
end

if(GroundTruthGiven)
	OutputS.GroundSeg = ResampledGroundCropped > 0;
	OutputS.GroundSeg = bwfill(OutputS.GroundSeg, 'holes', 8);
	OutputS.InitialDice = dices_coefficient(OutputS.InitialSeg, OutputS.GroundSeg);
	OutputS.FinalDice = dices_coefficient(OutputS.FinalSeg, OutputS.GroundSeg);
	
	FinalSegBoundaries = bwboundaries(OutputS.FinalSeg, 8);
	GroundSegBoundaries = bwboundaries(OutputS.GroundSeg, 8);
%	keyboard;
	if(length(FinalSegBoundaries) > 1 || length(GroundSegBoundaries) > 1)
		error('There should only be one boundary in the final and ground segmentations');
	end
	
% 	DX = bsxfun(@minus, FinalSegBoundaries{1}(:, 2), GroundSegBoundaries{1}(:, 2)');
% 	DY = bsxfun(@minus, FinalSegBoundaries{1}(:, 1), GroundSegBoundaries{1}(:, 1)');
% 	D = sqrt(DX .* DX + DY .* DY);
% 	OutputS.FinalHaussdorf = max(max(min(D)), max(min(D, [], 2)));
	OutputS.FinalHaussdorf = hausdorff_distance(OutputS.FinalSeg, OutputS.GroundSeg);
	
	%keyboard;
% 	clear DX DY D;
	
	if(DoThicknesses)
		if(~ThicknessCalculationErrorGround)
			[OutputS.GNDThickness.Contours, OutputS.GNDThickness.SolvedImage, OutputS.GNDThickness.InitialPoints, OutputS.GNDThickness.ArcLengthsLeft, OutputS.GNDThickness.ArcLengthsRight, OutputS.GNDThickness.LeftOffsetIDXStart, OutputS.GNDThickness.CurLeftOffsetIDX, OutputS.GNDThickness.RightOffsetIDXStart, OutputS.GNDThickness.CurRightOffsetIDX] = ...
			deal(GND.Contours, GND.SolvedImage, GND.InitialPoints, GND.ArcLengthsLeft, GND.ArcLengthsRight, GND.LeftOffsetIDXStart, GND.CurLeftOffsetIDX, GND.RightOffsetIDXStart, GND.CurRightOffsetIDX);

			[OutputS.GNDThickness.MaskInnerBoundary, ...
			OutputS.GNDThickness.MaskOuterBoundary, ...
			OutputS.GNDThickness.MaskFree, ...
			OutputS.GNDThickness.SolvedImage, ...
			OutputS.GNDThickness.XY, ...
			OutputS.GNDThickness.X, ...
			OutputS.GNDThickness.Y, ...
			OutputS.GNDThickness.NormFX, ...
			OutputS.GNDThickness.NormFY, ...
			OutputS.GNDThickness.ValidStreamlines, ...
			OutputS.GNDThickness.StartV] = deal(GND.MaskInnerBoundary, ...
			GND.MaskOuterBoundary, ...
			GND.MaskFree, ...
			GND.SolvedImage, ...
			GND.XY, ...
			GND.X, ...
			GND.Y, ...
			GND.NormFX, ...
			GND.NormFY, ...
			GND.ValidStreamlines, ...
			GND.StartV);

			OutputS.GND.FinalThicknesses = GND.ArcLengths;
		end
	end
end
	%[~, OutputS.NumRegions] = bwlabel(OutputS.FinalSeg);

save(OutputMAT, '-struct', 'OutputS');

function [H] = hausdorff_distance(A, B)

BoundaryA = bwboundaries(A);
BoundaryB = bwboundaries(B);

if(~isempty(BoundaryA) && ~isempty(BoundaryB))
	BoundaryA = cat(1, BoundaryA{:});
	BoundaryB = cat(1, BoundaryB{:});
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
	I = sub2ind(size(BW), IY, IX);
	I = unique(I);
	T = false(size(BW));
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

	JoiningSegments = imdilate(T, double(F > HalfMaximum));
	
% 		subplot 131;
% 		imshow(OldBW, []);
% 		subplot 132;
% 		imshow(BW, []);
% 		subplot 133;
% 		imshow(F > HalfMaximum, []);
% 		keyboard;

	clear I IX IY n ArcLength EndY StartY EndX StartX F T HalfMaximum SQ Angle AngleWeighting R;
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
	TemplateOverlap] = do_lk_and_otsu(ResampledAVW, MatchedFilter, MatchedFilterProb, MatchedFilterFornixProb, ResampledGroundAVW, real_total_offset, DoLK)

if nargin < 7
	DoLK = true;
end
InitialResampledAVW = ResampledAVW(real_total_offset(1):real_total_offset(1) + size(MatchedFilter, 1) - 1, real_total_offset(2):real_total_offset(2) + size(MatchedFilter, 2) - 1);
% take out bottom and top 2%

% S = sort(MatchedFilter(:));
% T = sort(InitialResampledAVW(:));
% SLow = S(round(length(S) * 0.02));
% SHigh = S(round(length(S) * 0.98));
% TLow = T(round(length(T) * 0.02));
% THigh = T(round(length(T) * 0.98));
% S = max(S, SLow);
% S = min(S, SLow);

% MatchedFilterClipped = max(MatchedFilter, SLow);
% MatchedFilterClipped = min(MatchedFilterClipped, SHigh);
% InitialResampledAVWClipped = max(InitialResampledAVW, TLow);
% InitialResampledAVWClipped = min(InitialResampledAVWClipped, THigh);
MatchedFilterRemapped = immatchhist(MatchedFilter, InitialResampledAVW);
% ResampledAVWClipped = max(ResampledAVW, TLow);
% ResampledAVWClipped = min(ResampledAVWClipped, THigh);

if(DoLK)
	[LKParameters, LKCost] = lk_weighted_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, MatchedFilter ./ max(MatchedFilter(:)), [0 0 0 0 real_total_offset(2) real_total_offset(1)], 100);
	%LKParameters = lk_weighted_run_affine_inv_comp(ResampledAVWClipped, MatchedFilterRemapped, MatchedFilterClipped, [0 0 0 0 real_total_offset(2) real_total_offset(1)], 100);
	%%
	
	%%
%	keyboard;
	%%
% 	subplot 121;
% 	imshow(InitialResampledAVW, []); colorbar;
% 	subplot 122;
% 	imshow(MatchedFilterRemapped, []); colorbar;
% 	
% 	%%
% 	LKParameters = lk_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, [0 0 0 0 real_total_offset(2) real_total_offset(1)], 100);
else
	LKParameters = [0 0 0 0 real_total_offset(2) real_total_offset(1)]';
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
%%
% TemplateLKIMGCC = bwconncomp(TemplateLKIMG > 0);
% TemplateLKIMGL = labelmatrix(TemplateLKIMGCC);
% R = regionprops(TemplateLKIMGCC, 'BoundingBox');
% if(length(R) > 1)
% 	disp('Bad image, returning');
% 	keyboard;
% 	return;
% end
ResampledAVWCropped = imcrop(ResampledAVW, TemplateLKIMGBoundingBox);
if(~isempty(ResampledGroundAVW))
	ResampledGroundCropped = imcrop(ResampledGroundAVW, TemplateLKIMGBoundingBox);
else
	ResampledGroundCropped = [];
end
TemplateProbLKIMGCropped = imcrop(TemplateProbLKIMG, TemplateLKIMGBoundingBox);
FornixProbLKIMGCropped = imcrop(FornixProbLKIMG, TemplateLKIMGBoundingBox);
%SIMGCropped = imcrop(SIMG, R.BoundingBox);
%AIMGCropped = imcrop(AIMG, R.BoundingBox);

[THRESH, SEG] = robust_otsu2(ResampledAVWCropped, [0.05, 0.98]);
SEG(ResampledAVWCropped > THRESH(end)) = 3;
%ResampledAVWCroppedForOtsu = (ResampledAVWCropped - min(ResampledAVWCropped(:))) ./ (max(ResampledAVWCropped(:)) - min(ResampledAVWCropped(:)));
%ResampledAVWCroppedForOtsu = round(ResampledAVWCroppedForOtsu * 255);
%T = otsu2_c(uint8(ResampledAVWCroppedForOtsu));
%T = double(T);
%OtsuMask = (ResampledAVWCroppedForOtsu >= T(end));
OtsuMask = SEG == 3;
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