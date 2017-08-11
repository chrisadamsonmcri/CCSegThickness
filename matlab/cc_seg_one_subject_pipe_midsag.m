function [ReturnCode] = cc_seg_one_subject_pipe_midsag(FileNameOrStruct, OutputDir, OutputPrefix, MidSagMethod)

% cc_seg_one_subject_pipe_seg(FileNameOrStruct, OutputDir, OutputPrefix, MidSagMethod)
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
%	OutputDir (char): directory to put output files
%	OutputPrefix (char): prefix of output files 
%	MidSagMethod (char): 'acpcdetect' or 'flirt'

% checks to see if we have a char array (so a nifti file) or a struct (CC image is a matlab array obtained in another function)
% what needs to come out of this
% AVWMidSag is the 2D midsagittal volume
% NIIPixdims [2]: X and Y dimensions of the image AVW
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%
% Changelog:
%	07/01/2014: included ACPCDETECT in with the m files, no need for users
%	to install separately

if(nargin < 5)
	MidSagMethod = 'acpcdetect';
end

% set up the export string for LD_LIBRARY_PATH that will be used to run FSL
% binaries
switch(computer)
	case 'GLNXA64'
		FSLLibraryString = 'export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; ';
	case 'GLNX86'
		FSLLibraryString = 'export LD_LIBRARY_PATH=/usr/lib/i386-linux-gnu:$LD_LIBRARY_PATH; ';
	case 'MACI64'
		FSLLibraryString = '';
	otherwise
		warning('Unsupported platform for FSL');
end

% set up the path for the acpcdetect binary
switch(computer)
	case {'GLNXA64', 'GLNX86'}
		ACPCDetectPath = 'bin';
	case 'MACI64'
		ACPCDetectPath = 'bin_maci';
end

% check for FSLDIR environment variable
FSLDIR = getenv('FSLDIR');
if(isempty(FSLDIR))
	warning('FSLDIR is not set');
end

% check for external binaries
%BinariesRequired = {'fslorient', 'fslswapdim'
[~, ~, ~] = mkdir(OutputDir);

%iptsetpref('ImShowAxesVisible', 'off');
ReturnCode = 0;
if(ischar(FileNameOrStruct))
	
	InFile = FileNameOrStruct;
	disp(['We are processing a file ' InFile]);
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
		
		if(false)%exist(OutputMAT, 'file') == 2 && strcmpi(MidSagMethod, 'acpcvdetect'))
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
			[status, result] = system([FSLLibraryString ' ' FSLDIR '/bin/fslorient -getorient ' InFile]);
			if(status ~= 0)
				disp('fslorient failed');
			end
			%disp(result);

			switch(deblank(result))
				case 'NEUROLOGICAL'
					OrientString = 'LR PA IS';
				case 'RADIOLOGICAL'
					OrientString = 'RL PA IS';
				otherwise
					error(['Unknown orientation retured by fslorient: ' deblank(result)]);
			end

			NIIFileForART = fullfile(ARTTempDirectory, 'in.nii');

			[~, ~] = system([FSLLibraryString ' ' FSLDIR '/bin/fslswapdim ' InFile ' ' OrientString ' ' NIIFileForART]);
			[NII, AVW] = load_nii(NIIFileForART);

			T = double(sort(AVW(:)));
			UU = floor(length(T) * 0.02):ceil(length(T) * 0.98);
			
			T = T(floor(length(T) * 0.02):ceil(length(T) * 0.98));
			%keyboard;
			OriginalRange = max(T(:)) - min(T(:));
			OriginalMin = min(T(:));
			
			T = (T - OriginalMin) ./ OriginalRange;
			T = uint8(round(T * 255));
			
			OtsuThresh = otsu1(T);
			OtsuThresh = OtsuThresh / 255 * OriginalRange + OriginalMin;
			clear T;
			%keyboard;
			OtsuSegOfAVW = false(size(AVW));
			OtsuSegOfAVW(AVW > OtsuThresh) = 1;
			%keyboard;
			%%
			%imshow(OtsuSegOfAVW(:, :, 30), []);
			
			%%
			OtsuSegOfAVWCC = bwconncomp(OtsuSegOfAVW, 26);
			L = labelmatrix(OtsuSegOfAVWCC);
			
			OtsuSegOfAVWAreas = cellfun('length', OtsuSegOfAVWCC.PixelIdxList);
			%keyboard;
			[~, I] = max(OtsuSegOfAVWAreas);
			
			[ID, JD, KD] = ind2sub(size(AVW), OtsuSegOfAVWCC.PixelIdxList{I});
			BoundingBoxI = [min(ID), max(ID)];
			BoundingBoxJ = [min(JD), max(JD)];
			BoundingBoxK = [min(KD), max(KD)];
			
			BoundingBoxK(1) = floor(BoundingBoxK(2) - 180 / NII.hdr.dime.pixdim(4));
			BoundingBoxK(1) = max(BoundingBoxK(1), 1);
			%keyboard;
			AVW = AVW(BoundingBoxI(1):BoundingBoxI(2), BoundingBoxJ(1):BoundingBoxJ(2), BoundingBoxK(1):BoundingBoxK(2));
			clear BoundingBoxI BoundingBoxJ BoundingBoxK ID JD KD I OtsuSegOfAVWAreas OtsuSegOfAVWCC;

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

			if(exist([NIIFileForART '.gz'], 'file') == 2)
				%gunzip([NIIFileForART '.gz']);
				delete([NIIFileForART '.gz']);
			end
			
			NIIFileARTOutput = fullfile(ARTTempDirectory, 'out_initial.nii');
			switch(lower(MidSagMethod))
				case 'acpcdetect'
					%keyboard;
					
					[MFilePath] = fileparts(which(mfilename));
					ARTHOME = fullfile(MFilePath, 'ART');
					
					CommandString = ['export ARTHOME=' ARTHOME '; ' fullfile(ARTHOME, ACPCDetectPath, 'acpcdetect') ' -i ' NIIFileForART ' -o ' NIIFileARTOutput];
					system(CommandString);
					%keyboard;
				case 'flirt'
					%NIIFileARTOutputAffine = fullfile(ARTTempDirectory, 'out_initial_affine.nii');
					NIIFileARTOutputAffineMat = fullfile(ARTTempDirectory, 'out_initial_affine.mat');
					D = 15;
					CommandString = [FSLLibraryString ' ' FSLDIR '/bin/flirt -in ' NIIFileForART ' -ref ' FSLDIR '/data/standard/MNI152_T1_1mm.nii.gz -dof 6 -searchrx ' num2str(-D) ' ' num2str(D) ' -searchry ' num2str(-D) ' ' num2str(D) ' -searchrz ' num2str(-D) ' ' num2str(D) ' -cost mutualinfo -out ' NIIFileARTOutput ' -omat ' NIIFileARTOutputAffineMat];
					system(CommandString);

				otherwise
					error('Midsagittal method not supported');
			end
			[NII, AVW] = load_nii(NIIFileARTOutput);
			AVWMidSag = extract_mid_slice(AVW);
			%keyboard;
			NIIPixdims = NII.hdr.dime.pixdim(3:4);
			save(OutputMAT, 'AVWMidSag', 'NIIPixdims');
			%NII.hdr.dime.dim(2) = 1;
			%NII.img = fliplr(AVWMidSag');
			%save_nii(NII, fullfile(OutputDir, [OutputPrefix '_midsag.nii.gz']));
			
			rmdir(ARTTempDirectory, 's');
            %ReturnCode = 0;
            %return;
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


OutputMAT = fullfile(OutputDir, [OutputPrefix '_midsag.mat']);

save(OutputMAT, 'AVWMidSag', 'NIIPixdims', 'MSPMethod');

function [AVWMidSag] = extract_mid_slice(AVW)

if(mod(size(AVW, 2), 2) == 0)
	disp('Even');
	XI = size(AVW, 2) / 2;
	XFrac = 0.5;
	AVWMidSag = AVW(size(AVW, 1):-1:1, XI, size(AVW, 3):-1:1) .* (1 - XFrac) + ...
	AVW(size(AVW, 1):-1:1, XI + 1, size(AVW, 3):-1:1) .* XFrac;
	clear XI XFrac MidsagSliceIDX;
else
	disp('Odd');
	XI = ceil(size(AVW, 2) / 2);
	AVWMidSag = AVW(size(AVW, 1):-1:1, XI, size(AVW, 3):-1:1);
	clear XI;
end
AVWMidSag = squeeze(AVWMidSag)';
AVWMidSag = double(AVWMidSag);
