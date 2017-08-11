function [ReturnCode] = cc_seg_one_subject_pipe_binmask(FileNameOrStruct, GroundTruthFile, OutputDir, OutputPrefix)

% cc_seg_one_subject_pipe_binmask(FileNameOrStruct, GroundTruthFile, OutputDir, OutputPrefix)
%
% DESCRIPTION
%	Performs preprocessing for the thickness calculation method given that we are getting a 2D CC mask as input.
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

% checks to see if we have a char array (so a nifti file) or a struct (CC image is a matlab array obtained in another function)
% what needs to come out of this
% AVWMidSag is the 2D midsagittal volume
% NIIPixdims [2]: X and Y dimensions of the image AVW

ReturnCode = 0;
if(ischar(FileNameOrStruct))
	%disp('We are processing a file');
	InFile = FileNameOrStruct;
	[NII, AVW] = load_nii(InFile);
	
	if(ismatrix(AVW))
		%disp('A 2D file');
		% 2D image
		AVWMidSag = double(AVW);
		AVWMidSag = padarray(AVWMidSag, [3, 3], 0, 'both');
		T = NII.hdr.dime.pixdim(2:4);
		NIIPixdims = T(NII.hdr.dime.dim(2:4) > 1);
		clear T;
	else
		% look for the only slice with
		error('Only 2D images supported');
	end
end
%return;

%if(exist(InFile, 'file') ~= 2)
%	error(['Input file ' InFile ' is not a regular file or does not exist']);
%end
%keyboard;

if(exist('cc_seg_atlases.mat', 'file') == 2)
	load('cc_seg_atlases.mat');
else
	AtlasDir = '/data/addo/cc_seg';
	[~, FornixProbAVW] = load_nii(fullfile(AtlasDir, 'all_fornix_prob.nii.gz'));
	[~, ProbAVW] = load_nii(fullfile(AtlasDir, 'all_cc_prob.nii.gz'));
	[TemplateNII, TemplateAVW] = load_nii(fullfile(AtlasDir, 'all_msp_mean.nii.gz'));
	save('cc_seg_atlases.mat', 'FornixProbAVW', 'ProbAVW', 'TemplateNII', 'TemplateAVW');
end
TemplatePixdims = TemplateNII.hdr.dime.pixdim(3:4);

AVWMidSagxx = (1:size(AVWMidSag, 2)) * NIIPixdims(1);
AVWMidSagxx = AVWMidSagxx - mean(AVWMidSagxx);

AVWMidSagyy = (1:size(AVWMidSag, 1)) * NIIPixdims(2);
AVWMidSagyy = AVWMidSagyy - mean(AVWMidSagyy);

[AVWMidSagX, AVWMidSagY] = meshgrid(AVWMidSagxx, AVWMidSagyy);

Templatexx = (1:size(TemplateAVW, 2)) * TemplatePixdims(1);
Templatexx = Templatexx - mean(Templatexx);

Templateyy = (1:size(TemplateAVW, 1)) * TemplatePixdims(2);
Templateyy = Templateyy - mean(Templateyy);

LeftLimit = max(min(Templatexx), min(AVWMidSagxx));
UpperLimit = max(min(Templateyy), min(AVWMidSagyy));
RightLimit = min(max(Templatexx), max(AVWMidSagxx));
LowerLimit = min(max(Templateyy), max(AVWMidSagyy));

[ResampleX, ResampleY] = meshgrid(LeftLimit:TemplatePixdims(1):RightLimit, UpperLimit:TemplatePixdims(2):LowerLimit);

ResampledAVW = interp2(AVWMidSagX, AVWMidSagY, AVWMidSag, ResampleX, ResampleY, 'nearest');
GivenSeg = ResampledAVW > 0;
GivenSeg = imdilate(GivenSeg, strel('disk', 2));

GivenSegCC = bwconncomp(GivenSeg);

if(GivenSegCC.NumObjects > 1)
	numPixels = cellfun(@numel,GivenSegCC.PixelIdxList);
	[~, idx] = max(numPixels);
	GivenSeg = false(size(GivenSeg));
	GivenSeg(GivenSegCC.PixelIdxList{idx}) = 1;
	GivenSegProps = regionprops(GivenSeg, 'BoundingBox');
	GivenSeg = imcrop(GivenSeg, GivenSegProps.BoundingBox);
end
[~, ~, ~] = mkdir(OutputDir);

OutputMAT = fullfile(OutputDir, [OutputPrefix '_seg.mat']);

OutputS.InitialSeg = GivenSeg;
OutputS.FinalSeg = GivenSeg;
OutputS.FinalSegArtefactsRemoved = GivenSeg;
OutputS.IMG = GivenSeg;
OutputS.TemplatePixdims = TemplatePixdims;
%keyboard;

save(OutputMAT, '-struct', 'OutputS');