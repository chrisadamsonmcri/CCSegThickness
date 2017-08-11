function [varargout] = cc_seg_process(InputDir, OutputDir, varargin)

% cc_seg_process(InputDir, OutputDir, param1, val1, param2, val2, ...)
%
% DESCRIPTION
%	CCSegPipeline processing script. Performs midsagittal alignment and
%	slice extraction, CC segmentation, thickness profile generation. The
%	script expected Gzipped NIFTI files (*.nii.gz) in the "InputDir" directory
%	and for each image creates mat files in the directory "OutputDir". Each
%	output file has a numeric prefix that is determined by the alphabetical
%	order of NIFTI files in "InputDir".
%
% PARAMETERS
%	InputDir (string): directory of compressed NIFTI files (*.nii.gz) for
%	processing
%	OutputDir (string): directory to write output files
%   REMAINING ARGUMENTS: execution steps
%		'midsag': perform midsagittal plane extraction
%		'seg': generate CC segmentations
%		'thickness': generate thickness profiles
%		'all': all processing steps, shorthand for 'midsag', 'seg', 'thickness'
%		'manedit_thickness': just generate thickness profiles of manedit
%		files, segmentations that have been manually edited
%	OPTIONAL PARAMETERS
%		The variables param1, val1, ... are parameter/value pairs. The
%		following parameter/value pairs are recognised:
%			'NumThicknessNodes' [1]: the number of nodes in the thickness
%			profile to generate, 100 if not given
%			'RedoDone' (logical) [1]: if false, do not perform analysis if
%			output files are already found. If true, always perform
%			analysis even if output file exists. false by default
%			'MSPMethod' (string): Midsagittal plane extraction method: 'acpcdetect'
%			(default) and 'flirt' are supported. 
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
%

% check arguments

if nargin < 3
	error('Not enough arguments');
end

if(exist(InputDir, 'dir') ~= 7)
	error('Input directory doesnt exist');
end

Modes = {'midsag', 'midsag_seg', 'seg', 'seg_thickness', 'midsag_seg_thickness', 'thickness', 'thickness_manedit'};

if(~ismember(Mode, Modes))
	error('Mode not supported');
end

DoMidSag = false;
DoSeg = false;
DoThickness = false;
DoThicknessManedit = false;

% 
% switch(lower(Mode))
% 	case 'seg'
% 		DoSeg = true;
% 	case 'seg_thickness'
% 		DoSeg = true;
% 		DoThickness = true;
% 	case 'seg_thickness'
% 		DoMidSag = true;
% 		DoSeg = true;
% 		DoThickness = true;
% 	case 'thickness'
% 		DoThickness = true;
% 	case 'thickness_manedit'
% 		DoThicknessManedit = true;
% end

% process options
NumThicknessNodes = 100;
RedoDone = false;
MSPMethod = 'acpcdetect';

z = 1;

while z <= length(varargin)
	switch(lower(varargin{z}))
		case {'numthicknessnodes', 'redodone', 'mspmethod'}
			if z == length(varargin)
				disp(['Argument ' varargin{z} ' has no associated value, ignoring']);
			else
				switch(lower(varargin{z}))
					case 'numthicknessnodes'
						NumThicknessNodes = varargin{z + 1};
					case 'redodone'
						RedoDone = varargin{z + 1};
					case 'mspmethod'
						MSPMethod = varargin{z + 1};
				end
			end
			z = z + 2;
		case {'midsag', 'seg', 'thickness', 'all', 'thickness_manedit'}
			switch(lower(varargin{z}))
				case 'midsag'
					DoMidSag = true;
				case 'seg'
					DoSeg = true;
				case 'thickness'
					DoThickness = true;
				case 'all'
					DoMidSag = true;
					DoSeg = true;
					DoThickness = true;
				case 'thickness_manedit'
					DoThicknessManedit = true;
			end
			z = z + 1;
		otherwise
			disp(['Argument ' varargin{z} ' unknown, ignoring']);
	end
end

for z = 1:2:length(varargin)
	if(length(varargin) >= z + 1)
		switch(lower(varargin{z}))
			case 'numthicknessnodes'
				NumThicknessNodes = varargin{z + 1};
			case 'redodone'
				RedoDone = varargin{z + 1};
			case 'mspmethod'
				MSPMethod = varargin{z + 1};
		end
	end
end

if(~isnumeric(NumThicknessNodes) || ~isscalar(NumThicknessNodes))
	error('NumThicknessNodes must be a numeric scalar');
end
if(~islogical(RedoDone) || ~isscalar(RedoDone))
	error('RedoDone must be a logical scalar');
end
if(~ischar(MSPMethod))
	error('MSPMethod must be a string');
end

MSPMethods = {'flirt', 'acpcdetect'};
if(~ismember(MSPMethod, MSPMethods))
	error(['Invalid MSPMethod, must be one of: ' deblank(sprintf('%s ', MSPMethods{:}))]);
end

if(NumThicknessNodes < 1 || floor(NumThicknessNodes) ~= NumThicknessNodes)
	error('NumThicknessNodes must be a positive integer');
end

D = dir(fullfile(InputDir, '*.nii.gz'));
[status, ~, ~] = mkdir(OutputDir);

if(status == false)
	error('Could not make output directory');
end
ExistingOutputFiles = findallfiles(OutputDir);
ExistingOutputFilePaths = cell(1, length(ExistingOutputFiles));
ExistingOutputFileNames = cell(1, length(ExistingOutputFiles));
ExistingOutputFileExtensions = cell(1, length(ExistingOutputFiles));

for z = 1:length(ExistingOutputFiles)
	[ExistingOutputFilePaths{z}, ExistingOutputFileNames{z}, ExistingOutputFileExtensions{z}] = fileparts(ExistingOutputFiles{z});
end

tokens = regexp(ExistingOutputFileNames, '^(\d\d\d\d)-(.*)$', 'tokens');

ExistingOutputFileNumbers = zeros(1, length(ExistingOutputFiles));
ExistingOutputFileSuffixes = cell(1, length(ExistingOutputFiles));

for z = 1:length(ExistingOutputFiles)
	ExistingOutputFileNumbers(z) = str2double(tokens{z}{1}{1});
	ExistingOutputFileSuffixes{z} = tokens{z}{1}{2};
end

IDX = 1:length(D);

[OutputBases{1:length(D)}] = deal(D(IDX).name);
OutputBases = strrep(OutputBases, '.nii.gz', '');

% for each existing output file, check to see where it is in the current
% list
ExistingOutputFileNewNumbers = zeros(1, length(ExistingOutputFiles));

for z = 1:length(OutputBases)
	T = strncmp(OutputBases{z}, ExistingOutputFileSuffixes, length(OutputBases{z}));
	F = find(T);
	if(~isempty(F))
		ExistingOutputFileNewNumbers(F) = z;
	end
	clear T F;
end

ExistingOutputFileNumbersNotEqual = find(ExistingOutputFileNumbers ~= ExistingOutputFileNewNumbers);

if(~isempty(ExistingOutputFileNumbersNotEqual))
	for z = 1:length(ExistingOutputFileNumbersNotEqual)
		OldFile = ExistingOutputFiles{ExistingOutputFileNumbersNotEqual(z)};
		NewFile = fullfile(...
			ExistingOutputFilePaths{ExistingOutputFileNumbersNotEqual(z)}, ...
			[sprintf('%04d', ExistingOutputFileNewNumbers(ExistingOutputFileNumbersNotEqual(z))), ...
			'-', ...
			ExistingOutputFileSuffixes{ExistingOutputFileNumbersNotEqual(z)}, ...
			ExistingOutputFileExtensions{ExistingOutputFileNumbersNotEqual(z)}]);
% 		disp([ExistingOutputFiles{ExistingOutputFileNumbersNotEqual(z)} ' -> ' ...
% 			fullfile(ExistingOutputFilePaths{ExistingOutputFileNumbersNotEqual(z)}, ...
% 			[sprintf('%04d', ExistingOutputFileNewNumbers(ExistingOutputFileNumbersNotEqual(z))) '-' ExistingOutputFileSuffixes{ExistingOutputFileNumbersNotEqual(z)} ExistingOutputFileExtensions{ExistingOutputFileNumbersNotEqual(z)}])]);
		movefile(OldFile, NewFile);
		clear OldFile NewFile;
	end
end

for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputPrefix = [sprintf('%04d', z) '-' BaseName];
	OutputBase = fullfile(OutputDir, OutputPrefix);
	
	if(DoMidSag)
		if((~RedoDone && exist([OutputBase '_midsag.mat'], 'file') ~= 2) || RedoDone)
            disp(OutputBase);
			cc_seg_one_subject_pipe_midsag(fullfile(InputDir, D(z).name), ...
				OutputDir, ...
				OutputPrefix, MSPMethod);
		end
	end
	
	if(DoSeg)
		if((~RedoDone && exist([OutputBase '_seg.mat'], 'file') ~= 2) || RedoDone)
            disp(OutputBase);
			cc_seg_one_subject_pipe_seg(OutputBase, []);
		end
	end
	if(DoThickness || DoThicknessManedit)
		if((~RedoDone && exist([OutputBase '_thickness.mat'], 'file') ~= 2) || RedoDone)
			if((DoThicknessManedit && exist([OutputBase '_seg_manedit.mat'], 'file') == 2) || DoThickness)
				cc_seg_one_subject_pipe_thickness(OutputBase, NumThicknessNodes);
			end
		end
	end
end
