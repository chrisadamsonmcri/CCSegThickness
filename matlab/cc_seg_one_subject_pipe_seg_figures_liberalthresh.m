function [ReturnCode] = cc_seg_one_subject_pipe_seg_figures(FileNameOrStruct, GroundTruthFile, OutputDir, OutputPrefix, MidSagMethod, DoAngleFigure)

% cc_seg_one_subject_pipe_seg(FileNameOrStruct, GroundTruthFile, OutputDir, OutputPrefix, MidSagMethod)
%
% DESCRIPTION
%	Performs corpus callosum segmentation the using illusory contours level set method and various mathematical morphology operations.
%	Performs the segmentation part only.
%
% PARAMETERS
%	First argument
%	NIIFileName (char): expected to be a NIFTI file
%	Data (struct) with fields:
%		.AVW: midsagittal slice
%		.PixDims: pixel dimensions (X, Y) of midasgittal slice
%	GroundTruthFile (char): NIFTI file that contains a ground truth
%	segmentation, use [] if no ground truth available
%	OutputDir (char): directory to put output files
%	OutputPrefix (char): prefix of output files 
%	MidSagMethod (char): 'acpcdetect' or 'flirt'

% checks to see if we have a char array (so a nifti file) or a struct (CC image is a matlab array obtained in another function)
% what needs to come out of this
% AVWMidSag is the 2D midsagittal volume
% NIIPixdims [2]: X and Y dimensions of the image AVW

if(nargin < 5)
	MidSagMethod = 'acpcdetect';
end

iptsetpref('ImShowAxesVisible', 'off');
ReturnCode = 0;
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
		[tokens] = regexp([name ext], '^(.+)((\.nii)(\.gz)?)$', 'tokens');
		%keyboard;
		OutputMAT = fullfile(pathstr, [tokens{1}{1} '_midsag.mat']);
		clear ext;
		clear pathstr name;
		
		if(exist(OutputMAT, 'file') == 2 && strcmpi(MidSagMethod, 'acpcdetect'))
			disp('Midsag already run, loading from file');
			load(OutputMAT);
		else
			disp(['Midsag running, method: ' MidSagMethod]);
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
%			keyboard;
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
			%cat(1, BoundingBoxI, BoundingBoxJ, BoundingBoxK)
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
			switch(lower(MidSagMethod))
				case 'acpcdetect'
					
					ACPCDetect = '/usr/local/ART/bin/acpcdetect';
					CommandString = ['export ARTHOME=/usr/local/ART; ' ACPCDetect ' -i ' NIIFileForART ' -o ' NIIFileARTOutput];
					system(CommandString);
				case 'flirt'
					FSLDIR = getenv('FSLDIR');
					if(isempty(FSLDIR))
						error('FSLDIR is not set');
					end
					NIIFileARTOutputAffine = fullfile(ARTTempDirectory, 'out_initial_affine.nii');
					NIIFileARTOutputAffineMat = fullfile(ARTTempDirectory, 'out_initial_affine.mat');
					D = 15;
					CommandString = ['export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; ' FSLDIR '/bin/flirt -in ' NIIFileForART ' -ref ' FSLDIR '/data/standard/MNI152_T1_1mm.nii.gz -dof 6 -searchrx ' num2str(-D) ' ' num2str(D) ' -searchry ' num2str(-D) ' ' num2str(D) ' -searchrz ' num2str(-D) ' ' num2str(D) ' -cost mutualinfo -out ' NIIFileARTOutput ' -omat ' NIIFileARTOutputAffineMat];
					system(CommandString);
% 					M = dlmread(NIIFileARTOutputAffineMat);
% 					[rotmat, skew, scales, transl, angles] = fsl_decompose_aff(M);
% 					MRigid = eye(4);
% 					MRigid(1:3, 1:3) = rotmat * skew(1:3, 1:3);
% 					MRigid(1:3, 4) = transl;
% 					RigidMatFile = fullfile(ARTTempDirectory, 'out_initial_rigid.mat');
% 					%NIIFileARTOutputRigid = fullfile(ARTTempDirectory, 'out_initial_rigid.nii');
% 					dlmwrite(RigidMatFile, MRigid, '\t');
% 					CommandString = ['export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; ' FSLDIR '/bin/flirt -in ' NIIFileForART ' -ref ' FSLDIR '/data/standard/MNI152_T1_1mm.nii.gz -applyxfm -init ' RigidMatFile ' -out ' NIIFileARTOutput];
% 					system(CommandString);
% 					
					%keyboard;
				otherwise
					error('Midsagittal method not supported');
			end
				
			%disp(CommandString);
			
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
			NII.hdr.dime.dim(2) = 1;
			NII.img = fliplr(AVWMidSag');
			save_nii(NII, fullfile(OutputDir, [OutputPrefix '_midsag.nii.gz']));
			%delete(fullfile(ARTTempDirectory, 'out.nii'));
			rmdir(ARTTempDirectory, 's');
			ReturnCode = 0;
			return;
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

OutputMAT = fullfile(OutputDir, [OutputPrefix '_seg.mat']);

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

if(exist('cc_seg_atlases.mat', 'file') == 2)
	load('cc_seg_atlases.mat');
else
	AtlasDir = '/data/addo/cc_seg';
	[~, FornixProbAVW] = load_nii(fullfile(AtlasDir, 'all_fornix_prob.nii.gz'));
	[~, ProbAVW] = load_nii(fullfile(AtlasDir, 'all_cc_prob.nii.gz'));
	[TemplateNII, TemplateAVW] = load_nii(fullfile(AtlasDir, 'all_msp_mean.nii.gz'));
	save('cc_seg_atlases.mat', 'FornixProbAVW', 'ProbAVW', 'TemplateNII', 'TemplateAVW');
end
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

% Already done figures for the template data

%keyboard;
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

[I, J] = find(WMSeg);
clear J;
%
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

% T1
F = gaussian_filter_max_1d(2);
WMSegSmoothed = imfilter(double(WMSeg), F(:), 'same', 'conv', 'replicate');
WMSegSmoothed = imfilter(WMSegSmoothed, F(:)', 'same', 'conv', 'replicate');

cc2 = normxcorr2(MatchedFilterProb, WMSegSmoothed);
% T2 FLAIR
% [F, G] = gaussian_filter_max_1d(2);
% WMSegSmoothed = imfilter(double(WMSeg), F(:), 'same', 'conv', 'replicate');
% WMSegSmoothed = imfilter(WMSegSmoothed, F(:)', 'same', 'conv', 'replicate');

% ax = imfilter(MatchedFilter, G(:)', 'same', 'conv', 'replicate');
% ay = imfilter(MatchedFilter, G(:), 'same', 'conv', 'replicate');
% MatchedFilterGradMAG = sqrt(ax .* ax + ay .* ay);
% 
% ax = imfilter(ResampledAVW, G(:)', 'same', 'conv', 'replicate');
% ay = imfilter(ResampledAVW, G(:), 'same', 'conv', 'replicate');
% ResampledAVWGradMAG = sqrt(ax .* ax + ay .* ay);
% 
%cc2 = normxcorr2(MatchedFilterGradMAG, ResampledAVWGradMAG);
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
%disp(['Expected location is of CC rank: ' num2str(ExpectedCCRank)]);
clear I;
%keyboard;

%FirstOffsetValid2R = find(yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1) & xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2), 1, 'first');
OffsetValid2R = find(...
	yoffset2R >= 1 & yoffset2R + size(MatchedFilter, 1) - 1 <= size(ResampledAVW, 1) & ...
	xoffset2R >= 1 & xoffset2R + size(MatchedFilter, 2) - 1 <= size(ResampledAVW, 2));

% figure for xcorr
%cc2 = normxcorr2(MatchedFilterProb, WMSegSmoothed);

SR = 1;
SC = 2;
clf;

TitleProps = {'FontSize', 20};
RecProps = {'Curvature', [0, 0], 'EdgeColor', 'w', 'LineWidth', 2};

AX = zeros(1, 2);

% AX(1) = subplot(SR, SC, 1);
% 
% MissingSize = size(cc2) - size(MatchedFilterProb);
% PS = floor(MissingSize / 2);
% 
% %T = padarray(MatchedFilterProb(2:end, 2:end), [1, 1], max(MatchedFilterProb(:)), 'both');
% T = padarray(MatchedFilterProb, PS, 0, 'both');
% 
% if(mod(MissingSize(1), 2) == 1)
% 	T = padarray(T, [1, 0], 0, 'post');
% end
% 
% if(mod(MissingSize(2), 2) == 1)
% 	T = padarray(T, [0, 1], 0, 'post');
% end
% 
% imshow(T, []);
% hold on;
% rectangle('Position', [PS(2), PS(1), size(MatchedFilterProb, 2), size(MatchedFilterProb, 1)], RecProps{:}); 
% title([roman_label(1), ' {\itPCC}'], TitleProps{:});

%AX(2) = subplot(SR, SC, 2);
AX(1) = subplot(SR, SC, 1);
PS = floor((size(MatchedFilterProb) - 1) / 2);

T = padarray(WMSegSmoothed, PS, 0, 'both');
if(mod(size(MatchedFilterProb, 1), 2) == 0)
	T = padarray(T, [1, 0], 0, 'post');
end
if(mod(size(MatchedFilterProb, 2), 2) == 0)
	T = padarray(T, [0, 1], 0, 'post');
end
imshow(T, []);
rectangle('Position', [PS(2), PS(1), size(WMSegSmoothed, 2), size(WMSegSmoothed, 1)], RecProps{:});
title([roman_label(1), ' {\itWMSeg}'], TitleProps{:});

%AX(3) = subplot(SR, SC, [3 4]);
AX(2) = subplot(SR, SC, 2);
%[newmap, newimg, cmapx] = bluewhitered_image(256, cc2);
T = max(cc2, 0);

T = (T - min(T(:))) ./ (max(T(:)) - min(T(:)));

newmap = repmat(linspace(0, 0.75, 256)', 1, 3); %gray(256);
newmap = flipdim(newmap, 1);

T = repmat(T, [1, 1, 3]);
T = interp1(linspace(0, 1, size(newmap, 1)), newmap(:, 1), T);
cmapx = linspace(0, 1, 256);

imshow(T, []);
hold on;
title([roman_label(2), ' {\itNCC}'], TitleProps{:});

AXPos = get(AX(2), 'Position');
C = colormap;

MidHeight = AXPos(2) + AXPos(4) / 2;
H = 0.2;
%plot(cc2RRegMaxJ, cc2RRegMaxI, 'm*');
F = find(OffsetValid2R);
plot(xoffset2R(OffsetValid2R(F(1))) + size(MatchedFilter, 2), yoffset2R(OffsetValid2R(F(1))) + size(MatchedFilter, 1), 'wx', 'LineWidth', 2, 'MarkerSize', 10);
plot(xoffset2R(OffsetValid2R(F(2))) + size(MatchedFilter, 2), yoffset2R(OffsetValid2R(F(2))) + size(MatchedFilter, 1), 'wx', 'LineWidth', 2, 'MarkerSize', 10);
plot(xoffset2R(OffsetValid2R(F(3:end))) + size(MatchedFilter, 2), yoffset2R(OffsetValid2R(F(3:end))) + size(MatchedFilter, 1), 'wo', 'LineWidth', 2);

plot(xoffset2RReg + size(MatchedFilter, 2), yoffset2RReg + size(MatchedFilter, 1), 'w+', 'MarkerSize', 10, 'LineWidth', 2);
axis on;
set(gca, 'XTick', [], 'YTick', []);

LeftAXPos = get(AX(1), 'Position');
RightAXPos = get(AX(2), 'Position');
RightAXPos(1) = LeftAXPos(1) + LeftAXPos(3) + 0.001;

set(AX(2), 'Position', RightAXPos);

AXPos = get(AX(2), 'Position');
LegAX = axes('Position', [AXPos(1) + AXPos(3) + 0.01, MidHeight - H / 2, AXPos(3) / 20, H]);
C = reshape(newmap, [size(C, 1), 1, 3]);
C = repmat(C, [1, 4]);
imagesc(linspace(0, (max(cmapx) - min(cmapx)) / (256 / 4)), cmapx, C);
axis xy;
set(gca, 'XTick', [], 'FontSize', 16, 'YAxisLocation', 'right');
ylabel('NCC', 'Rotation', 270, 'VerticalAlignment', 'bottom');
%keyboard;
[~, ~, ~] = mkdir(fullfile(OutputDir, 'paper_figures'));
OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-001-NCC.eps']);
%keyboard;
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off', 'Color', 'w');
% exportfig(gcf, OutputFile, 'Format', 'eps', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
cc_seg_save_figure_paper_eps(OutputFile);

% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG, 10);
% imwrite(IMG, OutputFile);
ReturnCode = 0;
%return;
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
	FornixProbLKIMG = cell(1, size(Offsets, 1));
	ResampledGroundCropped = cell(1, size(Offsets, 1));
	ResampledAVWCropped = cell(1, size(Offsets, 1));
	TemplateProbLKIMGCropped = cell(1, size(Offsets, 1));
	FornixProbLKIMGCropped = cell(1, size(Offsets, 1));
	OriginalOtsuMask = cell(1, size(Offsets, 1));
	OtsuMask = cell(1, size(Offsets, 1));
	TemplateOverlap = cell(1, size(Offsets, 1));
	LKOtsuSEG = cell(1, size(Offsets, 1));

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
	LKCCSeg = cell(1, length(TemplateOverlap));
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

		%for FIDX = 1:lenNumTries + 1gth(Props)
		%	MaskProps(z).(Props{FIDX}) = [TemplateProps.(Props{FIDX}) OtsuMaskProps.(Props{FIDX})];
		%end
		%clear TemplateProps OtsuMaskProps;
		T = TemplateProbLKIMGCropped{z} > 0.3;
		LKCCSeg{z} = OtsuMask{z} & ~FornixProbLKIMGCropped{z} > 0.3;
		LKCCSeg{z} = bwareaopen(LKCCSeg{z}, 50);
		DiceCoefficients(z) = dices_coefficient(T, LKCCSeg{z});
% 		HausdorffDistances(z) = hausdorff_distance(T, LKCCSeg{z});
% 		DeformationNorm(z) = norm(LKParameters{z}) + norm(LKParameters{z}([6, 5])' - Offsets(z, :));
		
		[Distances{z}, Angles{z}] = nearest_angle_distance(T, LKCCSeg{z}, DoAngleFigure && z == 1, fullfile(OutputDir, 'paper_figures', [OutputPrefix '-angles-distances.png']));
		SumDistances(z) = sum(Distances{z});
		SumAngles(z) = sum(abs(Angles{z}));

		TransformAngleX(z) = atan2(TY{z}(1, 2) - TY{z}(1, 1), TX{z}(1, 2) - TX{z}(1, 1));
		TransformAngleY(z) = atan2(TY{z}(2, 1) - TY{z}(1, 1), TX{z}(2, 1) - TX{z}(1, 1)) - pi / 2;
		clear T;
	end
	RigidFactor = abs(TransformAngleX - TransformAngleY);
	%keyboard;
	%%
	% clf;
	% for z = 1:length(OtsuMask)
	% 	subplot(2, 3, z);
	% 	imshow(LKCCSeg{z}, []);
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

	TestValues(1).text = 'LKCost';				TestValues(1).bestidx = MinLKCostI;				TestValues(1).values = LKCost;									TestValues(1).weight = 3;
	TestValues(2).text = 'CCTemplateCoverage';	TestValues(2).bestidx = MaxPCoverageI;			TestValues(2).values = PCoverage;								TestValues(2).weight = 1;
	TestValues(3).text = 'DiceScores';			TestValues(3).bestidx = MaxDiceCoefficientI;	TestValues(3).values = DiceCoefficients;						TestValues(3).weight = 1;
	%TestValues(3).text = 'Hausdorff';	TestValues(3).bestidx = MinHausdorffDistanceI;	TestValues(3).values = HausdorffDistances;
	TestValues(4).text = 'Distances';			TestValues(4).bestidx = MinDistancesI;			TestValues(4).values = SumDistances;							TestValues(4).weight = 1;
	TestValues(5).text = 'Angles';				TestValues(5).bestidx = MinAnglesI;				TestValues(5).values = SumAngles;								TestValues(5).weight = 1;
	TestValues(6).text = 'RigidFactor';			TestValues(6).bestidx = MinRigidFactorI;		TestValues(6).values = abs(TransformAngleX - TransformAngleY);	TestValues(6).weight = 1;
		
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
		SC = size(Offsets, 1);
		%keyboard;
		%clf;
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
			
			AX(2, z) = subplot(SR, SC, z + 3);
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
			'Parent', AX(2, 3));
	
		LegHandle = text(1.01, 0.5, {'Affine Alignment', '  Green - Initial', '  Red - After LKT'}, ...
 			'Units', 'Normalized', ...
 			'Color', 'k', ...
 			'HorizontalAlignment', 'left', ...
			'VerticalAlignment', 'middle', ...
 			'FontSize', 14, ...
			'Parent', AX(1, 3));
	
		for z = 1:size(AX, 2)
			N = 0.05;
			AXPos = get(AX(1, z), 'Position');
			AXPos(1) = AXPos(1) - N * z;
			set(AX(1, z), 'Position', AXPos);
			AXPos = get(AX(2, z), 'Position');
			AXPos(1) = AXPos(1) - N * z;
			set(AX(2, z), 'Position', AXPos);
		
		end
		
		if MainFittingIteration == 1
			OutputFile = fullfile(OutputDir, 'cc', [OutputPrefix '_cc_choice.png']);
		else
			OutputFile = fullfile(OutputDir, 'cc', [OutputPrefix '_cc_fallback.png']);
		end
		OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-002-LKCC.eps']);
		%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
		%keyboard;
		cc_seg_save_figure_paper_eps(OutputFile);
% 		FigPos = fullscreen_fig_pos;
% 		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% 		exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 		IMG = imread(OutputFile);
% 		IMG = imautocropwhite(IMG);
% 		imwrite(IMG, OutputFile);
	end
	%clf;
	%keyboard;
	% need to detect the situation where the LK transform has failed in all cases
	%return;
	if(~all(isinf(SumDistances)) && any(ValidIDX))
		%keyboard;
		%return;
		%%%% THE HAUSDORFF DISTANCE DOESNT ALWAYS WORK %%%%
		%%%% Try and get the rotation angles of the LKParameters, choose the one that rotates the least
		%BestAlignmentIDX = ValidIDXIDX(BestAlignmentIDX)
		InitialResampledAVW = InitialResampledAVW(ValidIDX);
		MatchedFilterRemapped = MatchedFilterRemapped(ValidIDX);
		LKParameters = LKParameters(ValidIDX);
		LKCCSeg = LKCCSeg(ValidIDX);
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
		LKCCSeg = LKCCSeg{BestLKI};
		TemplateOverlap = TemplateOverlap{BestLKI};
		break;
	else
		if MainFittingIteration == 2
			disp('Something went wrong, check your file');
			ReturnCode = 1;
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

OldLKCCSeg = LKCCSeg;
%keyboard;
LKCCSeg = LKCCSeg & (TemplateProbLKIMGCropped > 0.5);

clf;

SR = 1;
SC = 2;
AX = zeros(1, 2);
iptsetpref('ImShowAxesVisible', 'on');

TitleProps = {'FontSize', 20};
AX(1, 1) = subplot(SR, SC, 1); imshow(ResampledAVWCropped, []); show_segmentation(OldLKCCSeg); title([roman_label(1) ' {\itLKCCSeg*}'], TitleProps{:});
AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(LKCCSeg); title([roman_label(2) ' {\itCCSeg \leftarrow LKCCSeg* \cap PCC* > 0.5}'], TitleProps{:});
iptsetpref('ImShowAxesVisible', 'off');
set(AX, 'XTick', [], 'YTick', []);
OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-004-CCSegInit.png']);
%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
%keyboard;
cc_seg_save_figure_paper_png(OutputFile);
%FigPos = fullscreen_fig_pos;
%set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
%exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3) * 4, 'Height', FigPos(4) * 4, 'Color', 'rgb');
% 
%IMG = imread(OutputFile);
%IMG = imautocropwhite(IMG);
%imwrite(IMG, OutputFile);

%T = imerode(LKCCSeg, strel('disk', 1));
MeanIntensityInStartMask = mean(ResampledAVWCropped(LKCCSeg));
VarIntensityInStartMask = var(ResampledAVWCropped(LKCCSeg));

GaussianProb = normpdf(ResampledAVWCropped, MeanIntensityInStartMask, sqrt(VarIntensityInStartMask));
GaussianProb = GaussianProb ./ max(GaussianProb(:));
GaussianProb(ResampledAVWCropped > MeanIntensityInStartMask) = 1;

% LowPassMask = bwareaopen(OriginalOtsuMask, 200);
% T = ResampledAVWCropped;
% T(~LowPassMask) = NaN;
% ZC = colfilt(T, [15, 15], 'sliding', @mean_non_nan);
% 
% ResampledAVWCroppedLowPass = zeros(size(ResampledAVWCropped));
% ResampledAVWCroppedLowPass(LowPassMask) = ResampledAVWCropped(LowPassMask) ./ ZC(LowPassMask);
% T = imerode(LKCCSeg, strel('disk', 1));
% MeanIntensityInStartMaskLowPass = mean(ResampledAVWCroppedLowPass(T));
% VarIntensityInStartMaskLowPass = var(ResampledAVWCroppedLowPass(T));
% 
% GaussianProbLowPass = normpdf(ResampledAVWCroppedLowPass, MeanIntensityInStartMaskLowPass, sqrt(VarIntensityInStartMaskLowPass));
% GaussianProbLowPass = GaussianProbLowPass ./ max(GaussianProbLowPass(:));
% GaussianProb(ResampledAVWCropped > MeanIntensityInStartMask) = 1;
% clear T;
% 
% DoLowPassGraphics = false;
% if(DoLowPassGraphics)
% 	
% 	clf;
% 	SR = 4;
% 	SC = 2;
% 	subplot(SR, SC, 1);
% 	M = [min(ResampledAVWCropped(LKCCSeg)) max(ResampledAVWCropped(LKCCSeg))];
% 	imshow(ResampledAVWCropped, M);
% 	subplot(SR, SC, 2);
% 	M = [min(ResampledAVWCroppedLowPass(LKCCSeg)) max(ResampledAVWCroppedLowPass(LKCCSeg))];
% 	imshow(ResampledAVWCroppedLowPass, M);
% 	subplot(SR, SC, 3);
% 	RGB = zeros([size(GaussianProbLowPass), 3]);
% 	T = GaussianProb .* (ResampledAVWCropped < MeanIntensityInStartMask);
% 	RGB(:, :, 1) = T;
% 	T = GaussianProb .* (ResampledAVWCropped >= MeanIntensityInStartMask);
% 	RGB(:, :, 3) = T;
% 	imshow(RGB, []);
% 	imshow(GaussianProb, []);
% 	subplot(SR, SC, 4);
% 	RGB = zeros([size(GaussianProbLowPass), 3]);
% 	T = GaussianProbLowPass .* (ResampledAVWCroppedLowPass < MeanIntensityInStartMaskLowPass);
% 	RGB(:, :, 1) = T;
% 	T = GaussianProbLowPass .* (ResampledAVWCroppedLowPass >= MeanIntensityInStartMaskLowPass);
% 	RGB(:, :, 3) = T;
% 	%imshow(RGB, []);
% 	imshow(GaussianProbLowPass, []);
% 	subplot(SR, SC, 5);
% 	histfit(ResampledAVWCropped(LKCCSeg));
% 	subplot(SR, SC, 6);
% 	histfit(ResampledAVWCroppedLowPass(LKCCSeg));
% 	OutputFile = fullfile(OutputDir, 'lowpass', [OutputPrefix '.png']);
% 	
%  	[H, Bins] = hist(ResampledAVWCroppedLowPass(LKCCSeg), 100);
% 	F = gaussian_filter_max_1d(1);
% 	
% 	HD = imdilate(H(:), ones(7, 1));
% 	
%  	Peaks = (HD == H(:));
% % 	% find the highest peak
%  	PeaksIDX = find(Peaks);
%  	[PeakProb, I] = max(H(Peaks));
% 	MU = Bins(PeaksIDX(I));
% 	HF = conv(H(:), F(:), 'same');
% 	binwidth = Bins(2) - Bins(1);
% 	
% 	I = find(HF > PeakProb / 2, 1, 'last');
% 	IT = HF(I) - HF(I + 1);
% 	IFrac = (HF(I) - PeakProb / 2) / IT;
% 	FWHM = (Bins(I) + IFrac * binwidth - MU) * 2;
% 	SIGMA = FWHM ./ (2 * sqrt(2 * log(2)));
% 	%keyboard;
% % 	init = [Bins(PeaksIDX(I)) / 2, Bins(PeaksIDX(I))];
% % 	[label, model, llh] = emgm(ResampledAVWCroppedLowPass(LKCCSeg)', init);
%  	subplot(SR, SC, 7);
%  	bar(Bins, H);
%  	hold on;
%  	
%  	area = sum(LKCCSeg(:)) * binwidth;
% % 	y1 = normpdf(Bins, model.mu(1), model.Sigma(1)) * area;
% % 	y2 = normpdf(Bins, model.mu(2), model.Sigma(2)) * area;
% % 	plot(Bins, y1, Bins, y2);
% % 	%keyboard;
% 	%y1 = normpdf(Bins, MeanIntensityInStartMaskLowPass, sqrt(VarIntensityInStartMaskLowPass)) * area;
% 	y1 = normpdf(Bins, MU, SIGMA) * area;
% 	plot(Bins, y1);
% 	subplot(SR, SC, 8);
% 	imshow(normpdf(ResampledAVWCroppedLowPass, MU, SIGMA), []);
% 	[~, ~, ~] = mkdir(fullfile(OutputDir, 'lowpass'));
% 	%keyboard;
% 	FigPos = fullscreen_fig_pos;
% 	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% 	exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 	IMG = imread(OutputFile);
% 	IMG = imautocropwhite(IMG);
% 	imwrite(IMG, OutputFile);
% end
%return;
% ButterWorthLowerCutoff = 3;
% ButterWorthHigherCutoff = ButterWorthLowerCutoff + min(size(ResampledAVWCropped)) / 4;
% ResampledAVWCroppedButter = butterworthbpf(ResampledAVWCropped, ButterWorthLowerCutoff, ButterWorthHigherCutoff, 2);
% 
% MeanIntensityInStartMaskButter = mean(ResampledAVWCroppedButter(LKCCSeg));
% VarIntensityInStartMaskButter = var(ResampledAVWCroppedButter(LKCCSeg));
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
% [C, CC] = contour(LKCCSeg, [0.5, 0.5]);
% set(CC, 'Color', 'r');
% subplot(SR, SC, 2);
% imshow(ResampledAVWCroppedButter, []);
% subplot(SR, SC, 3);
% imshow(GaussianProb, []);
% subplot(SR, SC, 4);
% imshow(GaussianProbButter, []);
% subplot(SR, SC, 5);
% histfit(ResampledAVWCropped(LKCCSeg));
% subplot(SR, SC, 6);
% histfit(ResampledAVWCroppedButter(LKCCSeg));
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

SIGMA = 0.05;

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
%%
clf;

SR = 2;
SC = 2;
AX = zeros(2, 2);
iptsetpref('ImShowAxesVisible', 'on');

TitleProps = {'FontSize', 20};
AX(1, 1) = subplot(SR, SC, 1); imshow(ExclusionA, []); title([roman_label(1) ' {\itPenaltyA}'], TitleProps{:});
AX(1, 2) = subplot(SR, SC, 2); imshow(ExclusionB, []); title([roman_label(2) ' {\itPenaltyB}'], TitleProps{:});
AX(2, 1) = subplot(SR, SC, 3); imshow(ExclusionC, []); title([roman_label(3) ' {\itPenaltyC}'], TitleProps{:});
AX(2, 2) = subplot(SR, SC, 4); imshow(ExclusionWeight, []); title([roman_label(4) ' {\itPenaltyImage}'], TitleProps{:});

N = 0.075;
AXPos = get(AX(2, 1), 'Position'); AXPos(2) = AXPos(2) + N; set(AX(2, 1), 'Position', AXPos);
AXPos = get(AX(1, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; set(AX(1, 2), 'Position', AXPos);
AXPos = get(AX(2, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; AXPos(2) = AXPos(2) + N; set(AX(2, 2), 'Position', AXPos);
iptsetpref('ImShowAxesVisible', 'off');
set(AX, 'XTick', [], 'YTick', []);
OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-005-PenaltyImage']);
%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));

%%
%keyboard;
cc_seg_save_figure_paper_eps(OutputFile);
%cc_seg_save_figure_paper_png(OutputFile);
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);

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
% [C, CC] = contour(LKCCSeg, [0.5 0.5]);
% set(CC, 'Color', 'r');



%%%%%%%%%%%%%%%%%%%%%%%% Jump straight to  liberal threshold

%%%%%%%%%%%%%%%%%% Old method
% StartMaskBeforeSelect = LKCCSeg;
% [StrongI, StrongJ] = find(LKCCSeg);
% StartMask = bwselect(ExclusionWeight < 0.8, StrongJ, StrongI, 8);
% clear StrongI StrongJ;
% 
% %L = regionprops(StartMask, 'Area');
% %StartMaskAreas = [L.Area];
% StartMaskBeforeOpen = StartMask;
% StartMask = bwareaopen(StartMask, 200);
% 
% %subplot 122;
% % imshow(StartMask, []);
% %OriginalPSI = double(StartMask);
% %OriginalS = OtsuMask & (TemplateProbLKIMGCropped > 0.2) & (GaussianProb > 0.3 | ResampledAVWCropped > MeanIntensityInStartMask);
% %StartMask;
% 
% %T = S + FornixProbLKIMGCropped;
% %imshow(T, []);
% % subplot 321;
% % hold on;
% % [L, CC] = contour(S, [0.5, 0.5]);
% %keyboard;
% 
% %DMask = (GaussianProb > 0.2 | ResampledAVWCropped > MeanIntensityInStartMask) & ExclusionWeight < 0.5;
% %DMask = imdilate(DMask, strel('disk', 3));
% 
% % [L, CCT] = contour(DMask, [0.5, 0.5]);
% % set(CCT, 'Color', 'r');
% StartMaskBeforeFill = StartMask;
% 
% StartMask = bwfill(StartMask, 'holes', 8);
% %%
% clf;
% 
% SR = 2;
% SC = 2;
% AX = zeros(2, 2);
% iptsetpref('ImShowAxesVisible', 'on');
% 
% TitleProps = {'FontSize', 16};
% AX(1, 1) = subplot(SR, SC, 1); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeSelect); title([roman_label(1) ' {\itCCSeg}'], TitleProps{:});
% AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeOpen); title([roman_label(2) ' {\itCCSeg} \leftarrow SelectComponents({\itPenaltyImage} < 0.8, CCSeg)'], TitleProps{:});
% AX(2, 1) = subplot(SR, SC, 3); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeFill); title([roman_label(3) ' {\itCCSeg} \leftarrow RetainComponentsBySize({\itCCSeg}, 200)'], TitleProps{:});
% AX(2, 2) = subplot(SR, SC, 4); imshow(ResampledAVWCropped, []); show_segmentation(StartMask); title([roman_label(4) ' {\itCCSeg} \leftarrow FillHoles({\itCCSeg})'], TitleProps{:});
% 
% N = 0.05;
% AXPos = get(AX(2, 1), 'Position'); AXPos(2) = AXPos(2) + N; set(AX(2, 1), 'Position', AXPos);
% AXPos = get(AX(1, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; set(AX(1, 2), 'Position', AXPos);
% AXPos = get(AX(2, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; AXPos(2) = AXPos(2) + N; set(AX(2, 2), 'Position', AXPos);
% iptsetpref('ImShowAxesVisible', 'off');
% set(AX, 'XTick', [], 'YTick', []);
% OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-006-SelectOpenFill.png']);
% %[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
% cc_seg_save_figure_paper_png(OutputFile);
% %%
% %%
% % set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% % exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% % 
% % IMG = imread(OutputFile);
% % IMG = imautocropwhite(IMG);
% % imwrite(IMG, OutputFile);

StartMask = LKCCSeg;

StartMaskBeforeGauss = StartMask;
BW = ExclusionWeight < 1;

I = find(StartMask(:), 1, 'first');
[ID, JD] = ind2sub(size(StartMask), I);
BWS = bwselect(BW, JD, ID, 8);
AddedByReconstruct = BWS & ~StartMask;

StartMask = StartMask | AddedByReconstruct;

StartMaskBeforeClose = StartMask;

StartMask = imclose(StartMask, strel('disk', 4));

clf;

SR = 2;
SC = 2;
AX = zeros(2, 2);
iptsetpref('ImShowAxesVisible', 'on');

TitleProps = {'FontSize', 16};

AX(1, 1) = subplot(SR, SC, 1); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeGauss); title([roman_label(1) ' {\itCCSeg}'], TitleProps{:});
%AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeReconstruct); title({[roman_label(2) ' {\itCCSeg} \leftarrow'], 'FillHoles(LargestComponent(({\itCCSeg} \cap {\itNormProbCCSeg} > 0.01))'}, TitleProps{:});
AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeClose); title({[roman_label(2) ' {\itCCSeg} \leftarrow '] , ['SelectComponents({\itPenaltyImage} < 1 \cup {\itCCSeg}, {\itCCSeg})']}, TitleProps{:});
AX(2, 1) = subplot(SR, SC, 3); imshow(ResampledAVWCropped, []); show_segmentation(StartMask); title([roman_label(3) ' {\itCCSeg} \leftarrow {\itCCSeg} \bullet {\itD}_3'], TitleProps{:});
AX(2, 1) = subplot(SR, SC, 4); imshow(ResampledAVWCropped, []); show_segmentation(ExclusionWeight < 1); title([roman_label(4) ' {\itCCSeg} \leftarrow {\itCCSeg} \bullet {\itD}_3'], TitleProps{:});

% N = 0.05;
% AXPos = get(AX(2, 1), 'Position'); AXPos(2) = AXPos(2) + N; set(AX(2, 1), 'Position', AXPos);
% AXPos = get(AX(1, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; set(AX(1, 2), 'Position', AXPos);
% AXPos = get(AX(2, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; AXPos(2) = AXPos(2) + N; set(AX(2, 2), 'Position', AXPos);
iptsetpref('ImShowAxesVisible', 'off');
set(AX(AX ~= 0), 'XTick', [], 'YTick', []);

LeftAXPos = get(AX(1, 1), 'Position');
RightAXPos = get(AX(1, 2), 'Position');
BotAXPos = get(AX(2, 1), 'Position');

RightAXPos(1) = LeftAXPos(1) + LeftAXPos(3) + 0.035;
set(AX(1, 2), 'Position', RightAXPos);

BotAXPos(2) = LeftAXPos(2) - BotAXPos(4) - 0.05;
set(AX(2, 1), 'Position', BotAXPos);

keyboard;

clear StartMaskBefore*;


StartMaskCC = bwconncomp(StartMask, 8);
NumRegionsStartMask = StartMaskCC.NumObjects;
%[L, NumRegionsStartMask] = bwlabel(StartMask, 8);
%keyboard;
StartMaskBeforeFirstJoin = StartMask;
FirstJoinActivated = false;
if(NumRegionsStartMask > 1)
	FirstJoinActivated = true;
	T = bw_join_mst_by_boundaries(StartMask);
	StartMask = StartMask | (T & GaussianProb > 0.3);
	%DMask = DMask | (imdilate(T, strel('disk', 2)) & GaussianProb > 0.2);
	clear T;
end



% DMaskCC = bwconncomp(DMask);
% DMaskL = labelmatrix(DMaskCC);
% % 	imshow(label2rgb(DMaskL, 'lines', 'k'), []);
% % 	
% % find the region in DMask that mostly overlaps with the startmask, filter out other regions
% DMaskR = regionprops(DMaskCC, 'PixelIdxList');
% 
% Overlaps = zeros(1, length(DMaskR));
% I = find(StartMask);
% for z = 1:length(DMaskR)
% 	TF = ismember(DMaskR(z).PixelIdxList, I);
% 	Overlaps(z) = sum(TF) ./ length(I);
% 	clear TF;
% end
% clear I;
% 
% if(all(Overlaps < 0.975))
% 	clf;
% %	keyboard;
% 	%imshow(label2rgb(DMaskL, 'lines', 'k'), []);
% 	if(~any(DMask(:)))
% 		disp('Something went wrong, no mask to work with');
% 		ReturnCode = 1;
% 		return;
% 	end
% 	disp('overlaps dont match, the DMask needs to be joined');
% 	%NewDMask = bwareaopen(DMask, 200) & ismember(DMaskL, find(Overlaps > 0.1));
% 	%return;
% 	%NewDMaskJoiningSegments = bw_join_mst_by_boundaries(NewDMask);
% 	%NewDMaskJoined = NewDMask | NewDMaskJoiningSegments;
% 	%DMask = NewDMaskJoined;
% %	keyboard;
% 	try
% 		NewStartMaskJoiningSegments = bw_join_mst_by_boundaries(StartMask);
% 	catch e
% 		disp('acpcdetect failed');
% 		ReturnCode = 1;
% 		return;
% 	end
% 	disp('iajfoiadsfj');
% 	StartMask = StartMask | NewStartMaskJoiningSegments;
% 	
% 	%keyboard;
% else
% 	[~, L] = max(Overlaps);
% 	DMask = (DMaskL == L);
% end


StartMaskBeforeFill = StartMask;
StartMask = bwfill(StartMask, 'holes', 8);
StartMaskCC = bwconncomp(StartMask, 8);
NumRegionsStartMask = StartMaskCC.NumObjects;
StartMaskBeforeSecondJoin = StartMask;
SecondJoinActivated = false;
if(NumRegionsStartMask > 1)
	SecondJoinActivated = true;
	try
		NewStartMaskJoiningSegments = bw_join_mst_by_boundaries(StartMask);
	catch e
		disp('acpcdetect failed');
		ReturnCode = 1;
		return;
	end
	StartMask = StartMask | NewStartMaskJoiningSegments;
end


%keyboard;
clear Overlaps DMaskL L;

StartMaskBeforeClose = StartMask;
SE = strel('disk', 4);
%DMask = imclose(DMask, SE);
StartMask = imclose(StartMask, SE);

%%
clf;

SR = 2;
SC = 2;
AX = zeros(2, 2);
iptsetpref('ImShowAxesVisible', 'on');

TitleProps = {'FontSize', 16};
AX(1, 1) = subplot(SR, SC, 1); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeFirstJoin); title([roman_label(1) ' {\itCCSeg}'], TitleProps{:});
AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeSecondJoin); title({[roman_label(2) ' {\itCCSeg} \leftarrow {\itCCSeg} \cup '], 'JoiningSegments({\itCCSeg}) \cap {\itNormProbLKCCSeg} > 0.3'}, TitleProps{:});
AX(2, 1) = subplot(SR, SC, 3); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeClose); title([roman_label(3) ' {\itCCSeg} \leftarrow {\itCCSeg} \cup JoiningSegments({\itCCSeg})'], TitleProps{:});
AX(2, 2) = subplot(SR, SC, 4); imshow(ResampledAVWCropped, []); show_segmentation(StartMask); title([roman_label(4) ' {\itCCSeg} \leftarrow {\itCCSeg} \bullet {\itD}_4'], TitleProps{:});

N = 0.05;
AXPos = get(AX(2, 1), 'Position'); AXPos(2) = AXPos(2) + N; set(AX(2, 1), 'Position', AXPos);
AXPos = get(AX(1, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; set(AX(1, 2), 'Position', AXPos);
AXPos = get(AX(2, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; AXPos(2) = AXPos(2) + N; set(AX(2, 2), 'Position', AXPos);
iptsetpref('ImShowAxesVisible', 'off');
set(AX, 'XTick', [], 'YTick', []);
OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-007-MSTConnect.png']);
%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
cc_seg_save_figure_paper_png(OutputFile);
%%

% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);

clear StartMaskBefore*;
% 
% MeanIntensityInFinalStartMask = mean(ResampledAVWCropped(StartMask));
% VarIntensityInFinalStartMask = var(ResampledAVWCropped(StartMask));
% 
% GaussianProbFinal = normpdf(ResampledAVWCropped, MeanIntensityInFinalStartMask, sqrt(VarIntensityInFinalStartMask));
% GaussianProbFinal = GaussianProbFinal ./ max(GaussianProbFinal(:));
% %GaussianProbFinalNoPositive = GaussianProbFinal;
% GaussianProbFinal(ResampledAVWCropped > MeanIntensityInFinalStartMask) = 1;
% 
StartMaskBeforeGauss = StartMask;
% StartMask = StartMask & (GaussianProbFinal > 0.01);
% %%
% % subplot 221;
% % imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeGauss);
% % subplot 222;
% % imshow(ResampledAVWCropped, []); show_segmentation(GaussianProbFinal > 0.0003);
% % subplot 223;
% % imshow(GaussianProb, []);% show_segmentation(GaussianProbFinal > 0.0003);
% % subplot 224;
% % imshow(GaussianProbFinal, []);% show_segmentation(GaussianProbFinal > 0.0003);
% 
% %%
% 
% %DMask = DMask & (GaussianProbFinal > 0.001);
% 
% % DMaskCC = bwconncomp(DMask);
% % DMaskL = labelmatrix(DMaskCC);
% % DMaskAreas = cellfun('length', DMaskCC.PixelIdxList);
% % [~, L] = max(DMaskAreas);
% % %keyboard;
% % DMask = (DMaskL == L);
% 
% StartMaskBeforeFill = StartMask;
% StartMaskCC = bwconncomp(StartMask);
% StartMaskL = labelmatrix(StartMaskCC);
% StartMaskAreas = cellfun('length', StartMaskCC.PixelIdxList);
% [~, L] = max(StartMaskAreas);
% StartMask = (StartMaskL == L);
% 
% % DMask = bwfill(DMask, 'holes');
% StartMask = bwfill(StartMask, 'holes');

StartMaskBeforeReconstruct = StartMask;

DoReconstruct = true;

if(DoReconstruct)
% 		T = 5;
% 		M = min(T * double(StartMask), 6 - double(ExclusionWeight));
% 		G = imreconstruct(M, 6 - double(ExclusionWeight), 8);
% 		AddedByReconstruct = (G >= T) & ~StartMask;
% 		clear G M T;

	%	keyboard;

%	keyboard;
% FigPos = fullscreen_fig_pos;

%%

	%P = 0.5 * ExclusionB + 2 * ExclusionC;

% 	subplot 311; imshow(ResampledAVWCropped, []);
% 	subplot 312; imshow(P, []); show_segmentation(StartMaskBeforeSelect);
% 	subplot 313; imshow(P < 1.6, []);
	%FornixRegion = ExclusionA > 0.5;
	%SegFornixRegion = StartMask(FornixRegion);
	BW = ExclusionWeight < 1 | StartMask;
	I = find(StartMask(:), 1, 'first');
	[ID, JD] = ind2sub(size(StartMask), I);
	BWS = bwselect(BW, JD, ID, 8);
	AddedByReconstruct = BWS & ~StartMask;
	%clear I ID JD BWS BW;
% 	
% 	return;
%keyboard;
% 	%%
% 	subplot 121;
% 	imshow(AddedByReconstruct, []);
% 	subplot 122;
% 	imshow(AddedByBWSelect, []);
% 	
	
	%%
	StartMask = StartMask | AddedByReconstruct;
	%StartMask(FornixRegion) = SegFornixRegion;
	%StartMask = select_biggest_component(StartMask);
	
else
	AddedByReconstruct = false(size(StartMask));
end
%StartMask = bwfill(StartMask, 'holes');
StartMaskBeforeClose = StartMask;
StartMask = imclose(StartMask, strel('disk', 3));
%%
clf;

SR = 2;
SC = 2;
AX = zeros(2, 2);
iptsetpref('ImShowAxesVisible', 'on');

TitleProps = {'FontSize', 16};

AX(1, 1) = subplot(SR, SC, 1); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeGauss); title([roman_label(1) ' {\itCCSeg}'], TitleProps{:});
%AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeReconstruct); title({[roman_label(2) ' {\itCCSeg} \leftarrow'], 'FillHoles(LargestComponent(({\itCCSeg} \cap {\itNormProbCCSeg} > 0.01))'}, TitleProps{:});
AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeClose); title({[roman_label(2) ' {\itCCSeg} \leftarrow '] , ['SelectComponents({\itPenaltyImage} < 1 \cup {\itCCSeg}, {\itCCSeg})']}, TitleProps{:});
AX(2, 1) = subplot(SR, SC, 3); imshow(ResampledAVWCropped, []); show_segmentation(StartMask); title([roman_label(3) ' {\itCCSeg} \leftarrow {\itCCSeg} \bullet {\itD}_3'], TitleProps{:});

% N = 0.05;
% AXPos = get(AX(2, 1), 'Position'); AXPos(2) = AXPos(2) + N; set(AX(2, 1), 'Position', AXPos);
% AXPos = get(AX(1, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; set(AX(1, 2), 'Position', AXPos);
% AXPos = get(AX(2, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; AXPos(2) = AXPos(2) + N; set(AX(2, 2), 'Position', AXPos);
iptsetpref('ImShowAxesVisible', 'off');
set(AX(AX ~= 0), 'XTick', [], 'YTick', []);

LeftAXPos = get(AX(1, 1), 'Position');
RightAXPos = get(AX(1, 2), 'Position');
BotAXPos = get(AX(2, 1), 'Position');

RightAXPos(1) = LeftAXPos(1) + LeftAXPos(3) + 0.035;
set(AX(1, 2), 'Position', RightAXPos);

BotAXPos(2) = LeftAXPos(2) - BotAXPos(4) - 0.05;
set(AX(2, 1), 'Position', BotAXPos);

%%
%keyboard;
% 
% AX(1, 1) = subplot(SR, SC, 1); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeGauss); title([roman_label(1) ' {\itCCSeg}'], TitleProps{:});
% AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeReconstruct); title({[roman_label(2) ' {\itCCSeg} \leftarrow'], 'FillHoles(LargestComponent(({\itCCSeg} \cap {\itNormProbCCSeg} > 0.01))'}, TitleProps{:});
% AX(2, 1) = subplot(SR, SC, 3); imshow(ResampledAVWCropped, []); show_segmentation(StartMaskBeforeClose); title({[roman_label(3) ' {\itCCSeg} \leftarrow '] , ['SelectComponents({\itPenaltyImage} < 1 \cup {\itCCSeg}, {\itCCSeg})']}, TitleProps{:});
% AX(2, 2) = subplot(SR, SC, 4); imshow(ResampledAVWCropped, []); show_segmentation(StartMask); title([roman_label(4) ' {\itCCSeg} \leftarrow {\itCCSeg} \bullet {\itD}_3'], TitleProps{:});
% 
% N = 0.05;
% AXPos = get(AX(2, 1), 'Position'); AXPos(2) = AXPos(2) + N; set(AX(2, 1), 'Position', AXPos);
% AXPos = get(AX(1, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; set(AX(1, 2), 'Position', AXPos);
% AXPos = get(AX(2, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; AXPos(2) = AXPos(2) + N; set(AX(2, 2), 'Position', AXPos);
% iptsetpref('ImShowAxesVisible', 'off');
% set(AX, 'XTick', [], 'YTick', []);
% 

OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-008-GaussAndRecon.png']);
%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));

%%
cc_seg_save_figure_paper_png(OutputFile);
%keyboard;
%ReturnCode = 0;
%return;
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);

clear StartMaskBefore*;
%keyboard;

%% imreconstruct test
% subplot 121;
% imshow(StartMaskBeforeReconstruct, []);
% subplot 122;
% imshow(StartMask, []);
% imshow(6 - ExclusionWeight, []);
% subplot 222;
% imshow(G, []);
% subplot 223;
% imshow(AddedByReconstruct, []);
% keyboard;

%%
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
% OriginalPSI = zeros(size(StartMask));
% OriginalPSI(StartMask) = 1;
% OriginalPSI(~StartMask) = -1;
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
% NewPSI = OriginalPSI;
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
MorphGrad = T;
%BottomTopHat = BottomHat + TopHat;
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
clear T;
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
	OutlierPixels = false(size(OutlierPixelsStrong));
end

R = regionprops(double(FinalSeg), 'BoundingBox');
LeftColumn = floor(R.BoundingBox(1) + R.BoundingBox(3) * 2 / 3);
OutlierPixels(:, 1:LeftColumn) = 0;
clear R LeftColumn;

% dilate the outlier pixels
if(~any(OutlierPixels(:)))
	FinalSegOpenedM	= FinalSeg;
	VeinRemovalActivated = false;
else
	VeinRemovalActivated = true;
	SE = strel('disk', 3);
	
	% restrict the outlier pixels to the right third of the CC
	%keyboard;
	OutlierPixelsDilated = imdilate(OutlierPixels, SE);
	FinalSegWithout = FinalSeg & ~OutlierPixelsDilated;
	T = imdilate(OutlierPixelsDilated, ones(3, 3));
	FinalSegWithoutOutlierBorders = T & FinalSegWithout;
	% what you are doing here is
    % masking the part of the image that is just "in" the accepted area
    % surrounding the outliers, and then dilating them a bit
    % I guess this is 
	FinalSegOpened = imdilate(FinalSegWithoutOutlierBorders, SE) | FinalSegWithout;%imopen(FinalSegFirstCull, strel('disk', 3));
    %clear FinalSegWithoutOutlierBorders;
	FinalSegOpenedBeforeClosing = FinalSegOpened;
	% made the structuring element larger to fix it getting rid of the tail
	% of the splenium
	%FinalSegOpened = imclose(FinalSegOpened, strel('disk', 2));
	FinalSegOpened = imclose(FinalSegOpened, strel('disk', 4));
	%FinalSegOpenedBeforeMasking = FinalSegOpened;
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
	%%
	clf;

	SR = 3;
	SC = 3;
	AX = zeros(3, 3);
	iptsetpref('ImShowAxesVisible', 'on');

	TitleProps = {'FontSize', 12};
	AX(1, 1) = subplot(SR, SC, 1);
	imshow(ResampledAVWCropped, []);
	title([roman_label(1) ' {\itOriginalCCSeg}'], TitleProps{:});
	show_segmentation(FinalSeg);
% 	hold on;
% 	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
% 	set(CC, 'Color', 'r');
	
	AX(1, 2) = subplot(SR, SC, 2); imshow(MorphGrad, []); title([roman_label(2) ' {\itMorphGrad}'], TitleProps{:});
	
	AX(1, 3) = subplot(SR, SC, 3); imshow(ResampledAVWCropped, []); show_segmentation(OutlierPixels);
% 	[I, J] = find(OutlierPixels);
% 	hold on;
% 	plot(J, I, '*');
 	title([roman_label(3) ' {\itOutlierMask}'], TitleProps{:});
	
	AX(2, 1) = subplot(SR, SC, 4); imshow(ResampledAVWCropped, []); show_segmentation(FinalSegWithout); title({[roman_label(4) ' {\itCCSegNoOutliers} \leftarrow '], '{\itOriginalCCSeg} \cap', '({\itOutlierMask} \oplus {\itD}_{3})^{c} '}, TitleProps{:});
	AX(2, 2) = subplot(SR, SC, 5); imshow(ResampledAVWCropped, []); show_segmentation(FinalSegWithoutOutlierBorders); title({[roman_label(5) ' {\itCCSegNoOutliersAtBorders} \leftarrow '], '(({\itOutlierMask} \oplus {\itD}_{3}) \oplus {\itB}_{3}) \cap', '{\itCCSegNoOutliers}'}, TitleProps{:});
	AX(2, 3) = subplot(SR, SC, 6); imshow(ResampledAVWCropped, []); show_segmentation(FinalSegOpenedBeforeClosing); title({[roman_label(6) ' {\itCCSeg} \leftarrow '], '({\itCCSegNoOutliersAtBorders} \oplus {\itD}_3) \cup ', '{\itCCSegNoOutliers}'}, TitleProps{:});
	%AX(3, 2) = subplot(SR, SC, 8); imshow(FinalSegOpenedBeforeMasking, []); title([roman_label(8) ' {\itCCSeg} \leftarrow ({\itCCSegNoOutliersAtBorders} \oplus {\itD}_3) \cup {\itCCSegNoOutliers}'], TitleProps{:});
	AX(3, 2) = subplot(SR, SC, [7 8 9]); imshow(ResampledAVWCropped, []); show_segmentation(FinalSegOpenedM); title(roman_label(7), TitleProps{:});
	
	text(1.01, 0.5, {...
		'{\itCCSeg} \leftarrow {\itCCSeg} \bullet {\itD}_4', ...
		'{\itCCSeg} \leftarrow {\itCCSeg} \cap {\itOriginalCCSeg}', ...
		'{\itCCSeg} \leftarrow {\itCCSeg} \cup JoiningSegments({\itCCSeg})', ...
		'{\itCCSeg} \leftarrow LargestComponent({\itCCSeg})', ...
		'{\itCCSeg} \leftarrow {\itCCSeg} \bullet {\itD}_3', ...
		'{\itCCSeg} \leftarrow FillHoles({\itCCSeg})'}, ...
		'HorizontalAlignment', 'left', ...
		'VerticalAlignment', 'middle', ...
		'Units', 'normalized', ...
		'FontSize', 14);
		
	%({\itCCSegNoOutliersAtBorders} \oplus {\itD}_3) \cup {\itCCSegNoOutliers}'], TitleProps{:});
	%AX(2, 2) = subplot(SR, SC, 4); imshow(StartMask, []); title({[roman_label(4) ' {\itCCSeg} \leftarrow '] , ['SelectComponents({\itPenaltyImage} < 1 \cup {\itCCSeg}, {\itCCSeg})']}, TitleProps{:});

% 	N = 0.05;
% 	AXPos = get(AX(2, 1), 'Position'); AXPos(2) = AXPos(2) + N; set(AX(2, 1), 'Position', AXPos);
% 	AXPos = get(AX(1, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; set(AX(1, 2), 'Position', AXPos);
% 	AXPos = get(AX(2, 2), 'Position'); AXPos(1) = AXPos(1) - N * 1.25; AXPos(2) = AXPos(2) + N; set(AX(2, 2), 'Position', AXPos);
% 	iptsetpref('ImShowAxesVisible', 'off');
	
	
	for I = 1:2
		for J = 1:3
			N = 1.025;
			AXPos = get(AX(I, J), 'Position');
			AXPos(1) = AXPos(1) - (1 - N) * AXPos(3);
			AXPos(2) = AXPos(2) - (1 - N) * AXPos(4);
			AXPos(3:4) = AXPos(3:4) * N;
			set(AX(I, J), 'Position', AXPos);
		end
	end
		
	set(AX(AX ~= 0), 'XTick', [], 'YTick', []);
	TopAXPos = get(AX(2, 2), 'Position');
	AXPos = get(AX(3, 2), 'Position');
	AXPos([1 3]) = TopAXPos([1 3]);
	AXPos(2) = AXPos(2) + 0.05;
	set(AX(3, 2), 'Position', AXPos);
	OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-009-VeinRemoval.png']);
	%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
	cc_seg_save_figure_paper_png(OutputFile);
	%%
%	keyboard;
%  	FigPos = fullscreen_fig_pos;
%  	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
%  	exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3) * 4, 'Height', FigPos(4) * 4, 'Color', 'rgb');
% % 
%  	IMG = imread(OutputFile);
%  	IMG = imautocropwhite(IMG, 10);
%  	imwrite(IMG, OutputFile);

	clear StartMaskBefore*;
	
end

% FinalSegOpenedMOpened = imopen(FinalSegOpenedM, strel('disk', 1));
% FinalSegOpenedMOpenedLost = ~FinalSegOpenedMOpened & FinalSegOpenedM;
% 
% FinalSegOpenedMOpenedLostCC = bwconncomp(FinalSegOpenedMOpenedLost, 8);
% %FinalSegOpenedMOpenedLostL = labelmatrix(FinalSegOpenedMOpenedLostCC);
% FinalSegOpenedMThickened = FinalSegOpenedMOpened;
% for z = 1:length(FinalSegOpenedMOpenedLostCC.PixelIdxList)
% 	T = FinalSegOpenedM;
% 	T(FinalSegOpenedMOpenedLostCC.PixelIdxList{z}) = 0;
% 	TCC = bwconncomp(T);
% 	if(TCC.NumObjects > 1)
% 		T = false(size(FinalSegOpenedMOpened));
% 		T(FinalSegOpenedMOpenedLostCC.PixelIdxList{z}) = 1;
% 		FinalSegOpenedMThickened = FinalSegOpenedMThickened | T;%imdilate(T, strel('square', 3));
% 	end
% 	clear T TCC;
% end
% % %keyboard;
% 
% %%
% 
% %%
%clear FinalSegOpenedMOpened FinalSegOpenedMOpenedLostCC FinalSegOpenedMOpenedLost;
%OldFinalSegOpenedMThickened = FinalSegOpenedMThickened;
%OutlierPixels = bwareaopen(OutlierPixels, 6);
%OutlierPixels(:) = 0;%bwareaopen(OutlierPixels, 6);
% not doing the opening
FinalSegOpenedMThickened = FinalSegOpenedM;

%FinalSegOpenedMThickenedAfterOpening = FinalSegOpenedMThickened;
FinalSegOpenedMThickened = imdilate(FinalSegOpenedMThickened, strel('disk', 1));
FinalSegOpenedMThickened = bwfill(FinalSegOpenedMThickened, 'holes', 8);

FinalSeg = imdilate(FinalSeg, strel('disk', 1));
FinalSeg = bwfill(FinalSeg, 'holes', 8);

%%
clf;

SR = 1;
SC = 2;
AX = zeros(1, 2);
iptsetpref('ImShowAxesVisible', 'on');

TitleProps = {'FontSize', 14};
AX(1, 1) = subplot(SR, SC, 1);
imshow(ResampledAVWCropped, []); show_segmentation(FinalSegOpenedM);
title([roman_label(1) ' {\itCCSeg}'], TitleProps{:});
% hold on;
% [~, CC] = contour(FinalSeg, [0.5, 0.5]);
% set(CC, 'Color', 'r');

%AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(FinalSegOpenedMThickenedAfterOpening); title([roman_label(2) ' {\itCCSeg} \leftarrow {\itCCSeg} \o {\itB}_1'], TitleProps{:});
AX(1, 2) = subplot(SR, SC, 2); imshow(ResampledAVWCropped, []); show_segmentation(FinalSegOpenedMThickened); title([roman_label(2) ' {\itCCSeg} \leftarrow FillHoles({\itCCSeg} \oplus {\itD}_1)'], TitleProps{:});


% for I = 1:2
% 	for J = 1:3
% 		N = 1.025;
% 		AXPos = get(AX(I, J), 'Position');
% 		AXPos(1) = AXPos(1) - (1 - N) * AXPos(3);
% 		AXPos(2) = AXPos(2) - (1 - N) * AXPos(4);
% 		AXPos(3:4) = AXPos(3:4) * N;
% 		set(AX(I, J), 'Position', AXPos);
% 	end
% end
% 
set(AX(AX ~= 0), 'XTick', [], 'YTick', []);

LeftAXPos = get(AX(1), 'Position');
RightAXPos = get(AX(2), 'Position');
RightAXPos(1) = LeftAXPos(1) + LeftAXPos(3) + 0.001;

set(AX(2), 'Position', RightAXPos);

% TopAXPos = get(AX(2, 2), 'Position');
% AXPos = get(AX(3, 2), 'Position');
% AXPos([1 3]) = TopAXPos([1 3]);
% AXPos(2) = AXPos(2) + 0.05;
% set(AX(3, 2), 'Position', AXPos);
OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-010-FinalProcessing']);
%[~, ~, ~] = mkdir(fullfile(OutputDir, 'cc'));
cc_seg_save_figure_paper_png(OutputFile);
%%
% %keyboard;
% FigPos = fullscreen_fig_pos;
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);

%clear StartMaskBefore*;

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

PlotGraphics = false;

if(PlotGraphics)
	
	%%
	clf;
	SR = 3;
	SC = 3;
	
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
	
	[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	set(CC, 'Color', 'g');
	title('Original Image');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
	
	subplot(SR, SC, 5);
	imshow(ResampledAVWCropped, []);
	[I, J] = find(AddedByReconstruct);
	if(~isempty(I))
		hold on;
		plot(J, I, '*');
	end
	%	[~, CC] = contour(AddedByReconstruct, [0.5, 0.5]);
	%set(CC, 'Color', 'm');
	
	%[~, CC] = contour(OriginalPSI, [0, 0]);
	%set(CC, 'Color', 'b');
	%[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	%set(CC, 'Color', 'g');
	if(DoReconstruct)
		title('Original Image, added by reconstruct');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
	else
		title('Original Image, no reconstruct');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
	end

	subplot(SR, SC, 6);
	imshow(GaussianProb, []);
	hold on;
	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
	set(CC, 'Color', 'r');
	%[~, CC] = contour(OriginalPSI, [0, 0]);
	%set(CC, 'Color', 'b');
	[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	set(CC, 'Color', 'g');
	title('Gaussian probability of pixels');

	subplot(SR, SC, 7);
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

	subplot(SR, SC, 8);
	imshow(ExclusionWeight, []);
	hold on;
	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
	set(CC, 'Color', 'r');
	%[~, CC] = contour(OriginalPSI, [0, 0]);
	%set(CC, 'Color', 'b');
	[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
	set(CC, 'Color', 'g');
	title('Exclusion mask');

	if(GroundTruthGiven)
		subplot(SR, SC, 9);
		imshow(ResampledAVWCropped, []);	
		hold on;
		[~, CC] = contour(ResampledGroundCropped, [0.5, 0.5]);
		set(CC, 'Color', 'b');
		[~, CC] = contour(FinalSeg, [0.5, 0.5]);
		set(CC, 'Color', 'r');
		[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
		set(CC, 'Color', 'g');
		title({'Ground truth (blue)', 'Final (red)', 'Artefacts Removed (green)'});		

	end
	OutputPNG = fullfile(OutputDir, 'seg', [OutputPrefix '_seg.png']);
	%keyboard;
	[~, ~, ~] = mkdir(fullfile(OutputDir, 'seg'));
	
	FigPos = fullscreen_fig_pos;
	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
	exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
	
	IMG = imread(OutputPNG);
	IMG = imautocropwhite(IMG);
	imwrite(IMG, OutputPNG);
	%delete(gcf);
end

OutputS.InitialSeg = StartMask;
OutputS.FinalSeg = FinalSeg;
OutputS.FinalSegArtefactsRemoved = FinalSegOpenedMThickened;
OutputS.IMG = ResampledAVWCropped;
OutputS.TemplatePixdims = TemplatePixdims;
%keyboard;
OutputS.FirstJoinActivated = FirstJoinActivated;
OutputS.SecondJoinActivated = SecondJoinActivated;
OutputS.NumAddedByReconstruct = sum(AddedByReconstruct(:));
OutputS.VeinRemovalActivated = VeinRemovalActivated;

if(GroundTruthGiven)
	OutputS.GroundSeg = ResampledGroundCropped > 0;
	OutputS.GroundSeg = bwfill(OutputS.GroundSeg, 'holes', 8);
	OutputS.InitialDice = dices_coefficient(OutputS.InitialSeg, OutputS.GroundSeg);
	OutputS.FinalDice = dices_coefficient(OutputS.FinalSegArtefactsRemoved, OutputS.GroundSeg);

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
	OutputS.InitialHaussdorf = hausdorff_distance(OutputS.InitialSeg, OutputS.GroundSeg, TemplatePixdims);
	OutputS.FinalHaussdorf = hausdorff_distance(OutputS.FinalSegArtefactsRemoved, OutputS.GroundSeg, TemplatePixdims);
	
	%keyboard;
% 	clear DX DY D;

end
	%[~, OutputS.NumRegions] = bwlabel(OutputS.FinalSeg);

save(OutputMAT, '-struct', 'OutputS');

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

function [Distances, Angles, SmoothedBoundaryA, SmoothedBoundaryB] = nearest_angle_distance(A, B, DoFigure, OutputFigureFile)
% computes the sum of the differences in boundary vectors
% so for each boundary point in A, computes the angle between the tangent
% vector and the tangent vector of the nearest boundary point in B

%DoFigure = 1;

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
	SignT = sign(T);
	%TM = sum(TangentA .* -TangentB(DistanceIDX, :), 2);
	%TPLower = abs(TP)
	% clip
	%keyboard;
	T = abs(T);
	% clip to one since sometimes we have numerical precision problems
	T = min(T, 1);
	%T = max(T, -1);
	Angles = acos(T);
	if(DoFigure)
	%%
		clf;
		AX = axes;
		
%		SR = 1;%
		%SC = 3;
	% 	
	% 	AX = zeros(SR, SC);
	% 	AX(1) = subplot(SR, SC, 1);
	% 	imshow(A, []);
	% 	
	% 	AX(2) = subplot(SR, SC, 2);
	% 	imshow(B, []);
	% 	
	% 	AX(3) = subplot(SR, SC, 3);
		CroppingBox = [3.7195   24.9267  171.5122   69.1006];
		
		T = zeros(size(A));
		T(A & B) = 3;
		T(A & ~B) = 2;
		T(~A & B) = 1;
		T = imcrop(T, CroppingBox);
		imagesc(T);

		InsetCroppingBoxWidthHeight = [15.2378   6.3811];
		
		TopInsetCroppingBox = [102.8430, 7.4375, InsetCroppingBoxWidthHeight];
		LeftInsetCroppingBox = [50, 12, InsetCroppingBoxWidthHeight];
		
		FirstArrowLocation = [112, 9] + CroppingBox(1:2);
		SecondArrowLocation = [57, 16] + CroppingBox(1:2);
		
		% smoothedBoundaryA is the template
		XC = FirstArrowLocation(1) - SmoothedBoundaryA(:, 2);
		YC = FirstArrowLocation(2) - SmoothedBoundaryA(:, 1);
		[~, FirstArrowLocationIDX] = min(sqrt(XC .* XC + YC .* YC));
		FirstArrowLocationIDX = FirstArrowLocationIDX - 4;
		XC = SecondArrowLocation(1) - SmoothedBoundaryA(:, 2);
		YC = SecondArrowLocation(2) - SmoothedBoundaryA(:, 1);
		[~, SecondArrowLocationIDX] = min(sqrt(XC .* XC + YC .* YC));
		SecondArrowLocationIDX = SecondArrowLocationIDX - -2;
		IDXA = setdiff(1:size(SmoothedBoundaryA, 1), [FirstArrowLocationIDX, SecondArrowLocationIDX]);
		
		XC = SmoothedBoundaryA(FirstArrowLocationIDX, 2) - SmoothedBoundaryB(:, 2);
		YC = SmoothedBoundaryA(FirstArrowLocationIDX, 1) - SmoothedBoundaryB(:, 1);
		[~, FirstArrowPairIDX] = min(sqrt(XC .* XC + YC .* YC));
		
		XC = SmoothedBoundaryA(SecondArrowLocationIDX, 2) - SmoothedBoundaryB(:, 2);
		YC = SmoothedBoundaryA(SecondArrowLocationIDX, 1) - SmoothedBoundaryB(:, 1);
		[~, SecondArrowPairIDX] = min(sqrt(XC .* XC + YC .* YC));
		
		IDXB = setdiff(1:size(SmoothedBoundaryB, 1), [FirstArrowPairIDX, SecondArrowPairIDX]);
		QuiverScaling = 0.3;
		
		colormap gray;
		hold on;
		TAQH = quiver(SmoothedBoundaryA(IDXA, 2) - CroppingBox(1) + 0.5, SmoothedBoundaryA(IDXA, 1) - CroppingBox(2) + 0.5, TangentA(IDXA, 2), TangentA(IDXA, 1), QuiverScaling);
		set(TAQH, 'Color', 'r', 'LineWidth', 1.5);
		
		TBQH = quiver(SmoothedBoundaryB(IDXB, 2) - CroppingBox(1) + 0.5, SmoothedBoundaryB(IDXB, 1) - CroppingBox(2) + 0.5, TangentB(IDXB, 2), TangentB(IDXB, 1), QuiverScaling);
		set(TBQH, 'Color', 'b', 'LineWidth', 1.5);
		axis equal;
		set(gca, 'XLim', [1, size(T, 2)], 'YLim', [1, size(T, 1)], 'XTick', [], 'YTick', []);
		
		LegHandle = text(0.4, 0.25, {'Legend:', 'Dark Grey - {\itLKCCSeg}', 'Light Grey - {\itPCCMask}', 'White - {\itLKCCSeg} \cap {\itPCCMask}', '        BoundaryTangent({\itPCCMaskBP})', '        BoundaryTangent({\itLKCCSegBP})'}, ...
			'Units', 'Normalized', ...
			'Color', 'w', ...
			'HorizontalAlignment', 'left', ...
			'FontSize', 14);
		
		%QR = quiver(71.6, 58.3, 2, 0, 4);
		arrow([70.6, 58.3], [74.6, 58.3], 10, 'LineWidth', 2, 'EdgeColor', 'r', 'FaceColor', 'r', 'BaseAngle', 30);
		%set(QR, 'Color', 'r');
		arrow([70.6, 61.9], [74.6, 61.9], 10, 'LineWidth', 2, 'EdgeColor', 'b', 'FaceColor', 'b', 'BaseAngle', 30);
		%QR = quiver(70.6, 61.9, 2, 0, 4);
		%set(QR, 'Color', 'b');
		
		fix_figure_aspect(gcf, gca, size(T));
		ArrowProps = {20, 'BaseAngle', 30, 'Width', 2, 'EdgeColor', 'g', 'FaceColor', 'g'};
		DistanceArrowProps = {20, 'BaseAngle', 30, 'Width', 2, 'EdgeColor', 'y', 'FaceColor', 'y', 'Ends', 'both'};
		
		% do first arrow
		SX = SmoothedBoundaryA(FirstArrowLocationIDX, 2) - CroppingBox(1) + 0.5;
		SY = SmoothedBoundaryA(FirstArrowLocationIDX, 1) - CroppingBox(2) + 0.5;
		arrow([SX, SY], [SX + TangentA(FirstArrowLocationIDX, 2) * 4, SY + TangentA(FirstArrowLocationIDX, 1) * 4], ArrowProps{:});
		
		SX = SmoothedBoundaryB(FirstArrowPairIDX, 2) - CroppingBox(1) + 0.5;
		SY = SmoothedBoundaryB(FirstArrowPairIDX, 1) - CroppingBox(2) + 0.5;
		arrow([SX, SY], [SX + TangentB(FirstArrowPairIDX, 2) * 4, SY + TangentB(FirstArrowPairIDX, 1) * 4], ArrowProps{:});
		
		
		% Do second arrow
		SX = SmoothedBoundaryA(SecondArrowLocationIDX, 2) - CroppingBox(1) + 0.5;
		SY = SmoothedBoundaryA(SecondArrowLocationIDX, 1) - CroppingBox(2) + 0.5;
		
		arrow([SX, SY], [SX + TangentA(SecondArrowLocationIDX, 2) * 4, SY + TangentA(SecondArrowLocationIDX, 1) * 4], ArrowProps{:});
		
		TX = TangentB(SecondArrowPairIDX, 2) * 2;
		TY = TangentB(SecondArrowPairIDX, 1) * 2;
		if(SignT(SecondArrowLocationIDX) < 0)
			TX = -TX;
			TY = -TY;
		end
		
		SX = SmoothedBoundaryB(SecondArrowPairIDX, 2) - CroppingBox(1) + 0.5;
		SY = SmoothedBoundaryB(SecondArrowPairIDX, 1) - CroppingBox(2) + 0.5;
		arrow([SX, SY], [SX + TX * 2, SY + TY * 2], ArrowProps{:});
		
% 		TAQH = quiver(SmoothedBoundaryA(FirstArrowLocationIDX, 2) - CroppingBox(1) + 0.5, SmoothedBoundaryA(FirstArrowLocationIDX, 1) - CroppingBox(2) + 0.5, TangentA(FirstArrowLocationIDX, 2), TangentA(FirstArrowLocationIDX, 1), QuiverScaling * 5);
% 		set(TAQH, 'Color', 'y', 'LineWidth', 2);
% 		
		
		rectangle('Position', TopInsetCroppingBox, 'EdgeColor', 'm');
		rectangle('Position', LeftInsetCroppingBox, 'EdgeColor', 'm');
		
		InsetAX = zeros(2, 1);
		InsetAngleAX = zeros(2, 1);
		
		% do the top inset
		AXPos = get(AX, 'Position');
		InsetAX(1) = axes('Position', [AXPos(1) + AXPos(3) / 1.4, AXPos(2) + AXPos(4) / 1.1, AXPos(3) / 4, AXPos(4) / 4]);
		%annotation('rectangle', AXPos, 'Color', 'r');
		
		TTop = imcrop(T, TopInsetCroppingBox);
		imagesc(TTop);
		axis equal off;
		hold on;
		TAQH = quiver(SmoothedBoundaryA(IDXA, 2) - CroppingBox(1) - TopInsetCroppingBox(1) + 2, SmoothedBoundaryA(IDXA, 1) - CroppingBox(2) - TopInsetCroppingBox(2) + 2, TangentA(IDXA, 2), TangentA(IDXA, 1), QuiverScaling);
		P = {'MarkerSize', 3, 'MaxHeadSize', 2, 'LineWidth', 2};
		set(TAQH, 'Color', 'r', P{:});
		
		TBQH = quiver(SmoothedBoundaryB(IDXB, 2) - CroppingBox(1) - TopInsetCroppingBox(1) + 2, SmoothedBoundaryB(IDXB, 1) - CroppingBox(2) - TopInsetCroppingBox(2) + 2, TangentB(IDXB, 2), TangentB(IDXB, 1), QuiverScaling);
		set(TBQH, 'Color', 'b', P{:});
		%axis equal;
		xlim([1 size(TTop, 2)]);
		ylim([1 size(TTop, 1)]);
		%TAQH = quiver(SmoothedBoundaryA(FirstArrowLocationIDX, 2) - CroppingBox(1) - TopInsetCroppingBox(1) + 2, SmoothedBoundaryA(FirstArrowLocationIDX, 1) - CroppingBox(2) - TopInsetCroppingBox(2) + 2, TangentA(FirstArrowLocationIDX, 2), TangentA(FirstArrowLocationIDX, 1), QuiverScaling * 10);
		%set(TAQH, 'Color', 'y', 'LineWidth', 2, 'MarkerSize', 10, 'MaxHeadSize', 10);
		SX = SmoothedBoundaryA(FirstArrowLocationIDX, 2) - CroppingBox(1) - TopInsetCroppingBox(1) + 2;
		SY = SmoothedBoundaryA(FirstArrowLocationIDX, 1) - CroppingBox(2) - TopInsetCroppingBox(2) + 2;
		arrow([SX, SY], [SX + TangentA(FirstArrowLocationIDX, 2) * 2, SY + TangentA(FirstArrowLocationIDX, 1) * 2], ArrowProps{:});
		
		EX = SmoothedBoundaryB(FirstArrowPairIDX, 2) - CroppingBox(1) - TopInsetCroppingBox(1) + 2;
		EY = SmoothedBoundaryB(FirstArrowPairIDX, 1) - CroppingBox(2) - TopInsetCroppingBox(2) + 2;
		arrow([EX, EY], [EX + TangentB(FirstArrowPairIDX, 2) * 2, EY + TangentB(FirstArrowPairIDX, 1) * 2], ArrowProps{:});
		
		arrow([SX, SY], [EX, EY], DistanceArrowProps{:});
		
		text(7.1, 1.8, '{\itp}', 'Color', 'w', 'FontSize', 12);
		text(5.9, 6.3, '{\itq*}', 'Color', 'k', 'FontSize', 12);
		
		% do the tangents
		text(9.2, 2.1, '{\itT(p)}', 'Color', 'w', 'FontSize', 12);
		text(8.3, 5.7, '{\itT(q*)}', 'Color', 'w', 'FontSize', 12);
		
		text(6.1, 4.0, '$$||\overrightarrow{pq*}||_{2}$$', 'Interpreter', 'latex', 'Color', 'w', 'FontSize', 16, 'HorizontalAlignment', 'right');
		% draw dotted lines from the box in the image to the edges of the inset
		fix_figure_aspect(gcf, gca, size(TTop));
		AXPosMain = get(AX, 'Position');
		AXPosInset = get(InsetAX(1), 'Position');
		
		XLMain = xlim(AX);
		YLMain = ylim(AX);
		% find out the fraction of the figure width that the top left X point is
		XFracTopLeft = (TopInsetCroppingBox(1) - XLMain(1)) ./ (XLMain(2) - XLMain(1));
		YFracTopLeft = 1 - (TopInsetCroppingBox(2) - YLMain(1)) ./ (YLMain(2) - YLMain(1));
		
		XFracBotRight = (TopInsetCroppingBox(1) + TopInsetCroppingBox(3) - XLMain(1)) ./ (XLMain(2) - XLMain(1));
		YFracBotRight = 1 - (TopInsetCroppingBox(2) + TopInsetCroppingBox(4) - YLMain(1)) ./ (YLMain(2) - YLMain(1));
		
		annotation('line', [AXPosMain(1) + AXPosMain(3) * XFracTopLeft, AXPosInset(1)], [AXPosMain(2) + AXPosMain(4) * YFracTopLeft, AXPosInset(2) + AXPosInset(4)], 'Color', 'm', 'LineStyle', '--');
		annotation('line', [AXPosMain(1) + AXPosMain(3) * XFracBotRight, AXPosInset(1) + AXPosInset(3)], [AXPosMain(2) + AXPosMain(4) * YFracBotRight, AXPosInset(2)], 'Color', 'm', 'LineStyle', '--');
		
		AXPos = get(InsetAX(1), 'Position');
		InsetAngleAX(1) = axes('Position', [AXPos(1) + AXPos(3) / 2 - AXPos(3) / 4, AXPos(2) + AXPos(4) + 0.005, AXPos(3) / 2, AXPos(4) / 2]);
		QProps = {'MarkerSize', 10, 'LineWidth', 2};
		
		AngleArrowProps = {10, 'BaseAngle', 30, 'Width', 2, 'EdgeColor', 'k', 'FaceColor', 'k'};
		
		Q = quiver([0, 0], ...
			[0, 0], ...
			[TangentA(FirstArrowLocationIDX, 2) * 2, TangentB(FirstArrowPairIDX, 2) * 2], ...
			[TangentA(FirstArrowLocationIDX, 1) * 2, TangentB(FirstArrowPairIDX, 1) * 2], 2, QProps{:});
		
		axis(axis);
		arrow([0, 0], [TangentA(FirstArrowLocationIDX, 2) * 2, TangentA(FirstArrowLocationIDX, 1) * 2], AngleArrowProps{:});
		arrow([0, 0], [TangentB(FirstArrowPairIDX, 2) * 2, TangentB(FirstArrowPairIDX, 1) * 2], AngleArrowProps{:});
		arrow fixlimits;
		delete(Q);
		axis ij equal tight;
		
		%axis(axis);
		axis ij off equal;
		set(gca, 'XTick', [], 'YTick', []);
		%box off;
		text([1, 0.85, 2.03], [-0.22, 0.45, 0.1], {'{\itT(p)}', '{\itT(q*)}', '\angle{}{\itT(p)}{\itT(q*)}'}, 'Color', 'k', 'FontSize', 12);
		
		
		% do the left inset
		AXPos = get(AX, 'Position');
		InsetAX(2) = axes('Position', [AXPos(1) + AXPos(3) / 20, AXPos(2) + AXPos(4) / 1.1, AXPos(3) / 4, AXPos(4) / 4]);
		%annotation('rectangle', AXPos, 'Color', 'r');
		
		TLeft = imcrop(T, LeftInsetCroppingBox);
		imagesc(TLeft);
		axis equal off;
		hold on;
		TAQH = quiver(SmoothedBoundaryA(IDXA, 2) - CroppingBox(1) - LeftInsetCroppingBox(1) + 2, SmoothedBoundaryA(IDXA, 1) - CroppingBox(2) - LeftInsetCroppingBox(2) + 2, TangentA(IDXA, 2), TangentA(IDXA, 1), QuiverScaling);
		P = {'MarkerSize', 3, 'MaxHeadSize', 2, 'LineWidth', 2};
		set(TAQH, 'Color', 'r', P{:});
		
		TBQH = quiver(SmoothedBoundaryB(IDXB, 2) - CroppingBox(1) - LeftInsetCroppingBox(1) + 2, SmoothedBoundaryB(IDXB, 1) - CroppingBox(2) - LeftInsetCroppingBox(2) + 2, TangentB(IDXB, 2), TangentB(IDXB, 1), QuiverScaling);
		set(TBQH, 'Color', 'b', P{:});
		%axis equal;
		xlim([1 size(TLeft, 2)]);
		ylim([1 size(TLeft, 1)]);
		%TAQH = quiver(SmoothedBoundaryA(SecondArrowLocationIDX, 2) - CroppingBox(1) - LeftInsetCroppingBox(1) + 2, SmoothedBoundaryA(SecondArrowLocationIDX, 1) - CroppingBox(2) - LeftInsetCroppingBox(2) + 2, TangentA(SecondArrowLocationIDX, 2), TangentA(SecondArrowLocationIDX, 1), QuiverScaling * 10);
		%set(TAQH, 'Color', 'y', 'LineWidth', 2, 'MarkerSize', 10, 'MaxHeadSize', 10);
		SX = SmoothedBoundaryA(SecondArrowLocationIDX, 2) - CroppingBox(1) - LeftInsetCroppingBox(1) + 2;
		SY = SmoothedBoundaryA(SecondArrowLocationIDX, 1) - CroppingBox(2) - LeftInsetCroppingBox(2) + 2;
		arrow([SX, SY], [SX + TangentA(SecondArrowLocationIDX, 2) * 2, SY + TangentA(SecondArrowLocationIDX, 1) * 2], ArrowProps{:});
		
		EX = SmoothedBoundaryB(SecondArrowPairIDX, 2) - CroppingBox(1) - LeftInsetCroppingBox(1) + 2;
		EY = SmoothedBoundaryB(SecondArrowPairIDX, 1) - CroppingBox(2) - LeftInsetCroppingBox(2) + 2;
		TX = TangentB(SecondArrowPairIDX, 2) * 2;
		TY = TangentB(SecondArrowPairIDX, 1) * 2;
		if(SignT(SecondArrowLocationIDX) < 0)
			TX = -TX;
			TY = -TY;
		end
		arrow([EX, EY], [EX + TX, EY + TY], ArrowProps{:});
		arrow([SX, SY], [EX, EY], DistanceArrowProps{:});
		text(9.4, 3.7, '$$||\overrightarrow{pq*}||_{2}$$', 'Interpreter', 'latex', 'Color', 'w', 'FontSize', 16, 'HorizontalAlignment', 'left');
		
		text([11.4, 7.1], [5.9, 4.3], {'{\itp}', '{\itq*}'}, 'Color', 'w', 'FontSize', 12);
		% do the tangents
		text(13.2, 5.2, '{\itT(p)}', 'Color', 'w', 'FontSize', 12);
		text(8.3, 1.8, '{\itT(q*)}', 'Color', 'w', 'FontSize', 12);
		
		% draw dotted lines from the box in the image to the edges of the inset
		fix_figure_aspect(gcf, gca, size(TLeft));
		
		AXPosMain = get(AX, 'Position');
		AXPosInset = get(InsetAX(2), 'Position');
		
		XLMain = xlim(AX);
		YLMain = ylim(AX);
		% find out the fraction of the figure width that the Left left X point is
		XFracTopRight = (LeftInsetCroppingBox(1) + LeftInsetCroppingBox(3) - XLMain(1)) ./ (XLMain(2) - XLMain(1));
		YFracTopRight = 1 - (LeftInsetCroppingBox(2) - YLMain(1)) ./ (YLMain(2) - YLMain(1));
		
		XFracBotLeft = (LeftInsetCroppingBox(1) - XLMain(1)) ./ (XLMain(2) - XLMain(1));
		YFracBotLeft = 1 - (LeftInsetCroppingBox(2) + LeftInsetCroppingBox(4) - YLMain(1)) ./ (YLMain(2) - YLMain(1));
		
		annotation('line', [AXPosMain(1) + AXPosMain(3) * XFracTopRight, AXPosInset(1) + AXPosInset(3)], [AXPosMain(2) + AXPosMain(4) * YFracTopRight, AXPosInset(2) + AXPosInset(4)], 'Color', 'm', 'LineStyle', '--');
		annotation('line', [AXPosMain(1) + AXPosMain(3) * XFracBotLeft, AXPosInset(1)], [AXPosMain(2) + AXPosMain(4) * YFracBotLeft, AXPosInset(2)], 'Color', 'm', 'LineStyle', '--');
		
		AXPos = get(InsetAX(2), 'Position');
		InsetAngleAX(2) = axes('Position', [AXPos(1) + AXPos(3) / 2 - AXPos(3) / 8, AXPos(2) + AXPos(4) + 0.01, AXPos(3) / 4, AXPos(4) / 2]);
		
		Q = quiver([0, 0], ...
 			[0, 0], ...
 			[TangentA(SecondArrowLocationIDX, 2) * 4, TX * 2], ...
 			[TangentA(SecondArrowLocationIDX, 1) * 4, TY * 2], 0, QProps{:});
		axis(axis);
		arrow([0, 0], [TangentA(SecondArrowLocationIDX, 2) * 4, TangentA(SecondArrowLocationIDX, 1) * 4], AngleArrowProps{:});
		arrow([0, 0], [TX * 2, TY * 2], AngleArrowProps{:});
		arrow fixlimits;
		delete(Q);
		axis ij equal off tight;
		%set(gca, 'XTick', [], 'YTick', []);
		YL = ylim;
		YL(1) = YL(1) + 0.01;
		ylim(YL);
		%box off;
		text([2, -3, 0.5], [0, TY, -2], {'{\itT(p)}', '{\itT(q*)}', '\angle{}{\itT(p)}{\itT(q*)}'}, 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'middle')
		%text(1, -2, '\angle{}{\itp}{\itq*}', 'Color', 'k', 'FontSize', 12);
		
	%%
		[pathstr] = fileparts(OutputFigureFile);
		[~, ~, ~] = mkdir(pathstr);
		%OutputFile = fullfile(OutputDir, 'paper_figures', [OutputPrefix '-001-NCC.png']);
		%keyboard;
		%cc_seg_save_figure_paper_eps(OutputFigureFile);
		cc_seg_save_figure_paper_png(OutputFigureFile);
% 		FigPos = fullscreen_fig_pos;
% 		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off', 'Color', 'w');
% 		exportfig(gcf, OutputFigureFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 		IMG = imread(OutputFigureFile);
% 		IMG = imautocropwhite(IMG, 10);
% 		imwrite(IMG, OutputFigureFile);

%		keyboard;
	end
	
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
InitialResampledAVW = ResampledAVW(real_total_offset(1):real_total_offset(1) + size(MatchedFilter, 1) - 1, real_total_offset(2):real_total_offset(2) + size(MatchedFilter, 2) - 1);
% take out bottom and top 2%

% [F, G] = gaussian_filter_max_1d(1);
% ax = imfilter(ResampledAVW, G(:)', 'same', 'conv', 'replicate');
% ay = imfilter(ResampledAVW, G(:), 'same', 'conv', 'replicate');
% ResampledAVWGradMAG = sqrt(ax .* ax + ay .* ay);

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

% ax = imfilter(MatchedFilter, G(:)', 'same', 'conv', 'replicate');
% ay = imfilter(MatchedFilter, G(:), 'same', 'conv', 'replicate');
% MatchedFilterGradMAG = sqrt(ax .* ax + ay .* ay);

MatchedFilterRemapped = immatchhist(MatchedFilter, InitialResampledAVW);
%MatchedFilterGradMAGRemapped = immatchhist(MatchedFilterGradMAG, ResampledAVWGradMAG);
% ResampledAVWClipped = max(ResampledAVW, TLow);
% ResampledAVWClipped = min(ResampledAVWClipped, THigh);


if(DoLK)
	%keyboard;
	[LKParameters, LKCost] = lk_weighted_run_affine_inv_comp(ResampledAVW, MatchedFilterRemapped, MatchedFilter ./ max(MatchedFilter(:)), [0 0 0 0 real_total_offset(2) real_total_offset(1)], 100);
	%[LKParameters, LKCost] = lk_weighted_run_affine_inv_comp(ResampledAVWGradMAG, MatchedFilterGradMAGRemapped, MatchedFilter ./ max(MatchedFilter(:)), [0 0 0 0 real_total_offset(2) real_total_offset(1)], 100);	
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
	LKCost = 0;
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

[THRESH, OtsuSEG] = robust_otsu2(ResampledAVWCropped, [0.05, 0.98]);
OtsuSEG(ResampledAVWCropped > THRESH(end)) = 3;
%ResampledAVWCroppedForOtsu = (ResampledAVWCropped - min(ResampledAVWCropped(:))) ./ (max(ResampledAVWCropped(:)) - min(ResampledAVWCropped(:)));
%ResampledAVWCroppedForOtsu = round(ResampledAVWCroppedForOtsu * 255);
%T = otsu2_c(uint8(ResampledAVWCroppedForOtsu));
%T = double(T);
%OtsuMask = (ResampledAVWCroppedForOtsu >= T(end));
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

function [IMGH, CH] = show_segmentation(IMG)

hold on;
T = cat(3, double(IMG), zeros(size(IMG)), zeros(size(IMG)));
IMGH = imagesc(T, 'AlphaData', double(IMG) * 0.3);
[~, CH] = contour(double(IMG), [0.5, 0.5]);
set(CH, 'Color', 'r', 'LineWidth', 4);
%set(CH, 'Color', 'r', 'LineWidth', 1);

function [OBW] = select_biggest_component(BW)

CC = bwconncomp(BW);
CCAreas = cellfun('length', CC.PixelIdxList);
[~, I] = max(CCAreas);
OBW = false(size(BW));
%keyboard;
OBW(CC.PixelIdxList{I}) = 1;