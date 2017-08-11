function [ReturnCode] = cc_seg_mats_to_hdf5s(InputDir, OutputDir)

% converts all the mat files output by matlab as mat files in InputDir to
% the corresponding hdf5 files that would be output by the python version
% in OutputDir

if(exist(InputDir, 'dir') ~= 7)
	error('InputDir is not a directory');
end

if(exist(OutputDir, 'dir') ~= 7)
	disp('OutputDir is not a directory, making');
	[~, ~, ~] = mkdir(OutputDir);
end

D = dir(fullfile(InputDir, '*.mat'));

[FileNames{1:length(D)}] = deal(D.name);

tokens = regexp(FileNames, '^(.*)_(midsag|seg|seg_manedit|thickness)\.mat$', 'tokens');

for z = 1:length(FileNames)
	[pathstr, name, ext] = fileparts(FileNames{z});
	OutFile = fullfile(OutputDir, [name(6:end) '.hdf5']);
	%disp(OutFile);
	if(exist(OutFile, 'file') == 2)
		delete(OutFile);
	end
	switch(tokens{z}{1}{2})
		case 'midsag'
			T = load(fullfile(InputDir, FileNames{z}));
			if(isfield(T, 'NIIPixdims'))
				h5create(OutFile, '/NIIPixdims', size(T.NIIPixdims));
				h5write(OutFile, '/NIIPixdims', T.NIIPixdims);
			end
			if(isfield(T, 'AVWMidSag'))
				write_hdf_image(OutFile, '/midSagAVW', T.AVWMidSag);
				
			end
		case 'seg'
			T = load(fullfile(InputDir, FileNames{z}));
			%keyboard;
			
			if(isfield(T, 'FinalSeg'))
				write_hdf_image(OutFile, '/finalSeg', T.FinalSeg);
			end
			
			if(isfield(T, 'FinalSegArtefactsRemoved'))
				write_hdf_image(OutFile, '/finalSegNoArtefacts', T.FinalSegArtefactsRemoved);
			end
			
			if(isfield(T, 'InitialSeg'))
				write_hdf_image(OutFile, '/initialSeg', T.InitialSeg);
			end
			
			if(isfield(T, 'IMG'))
				write_hdf_image(OutFile, '/IMG', T.IMG);
			end
			
			if(isfield(T, 'TemplatePixdims'))
				h5create(OutFile, '/templatePixdims', size(T.TemplatePixdims));
				h5write(OutFile, '/templatePixdims', T.TemplatePixdims);
			end 
		case 'seg_manedit'
			T = load(fullfile(InputDir, FileNames{z}));
			if(isfield(T, 'FinalSegManEdit'))
				write_hdf_image(OutFile, '/finalSegManEdit', T.FinalSegManEdit);
			end
			
		case 'thickness'
			T = load(fullfile(InputDir, FileNames{z}));
	end
	clear T;
	%keyboard;
end

function [varargout] = write_hdf_image(OutFile, DatasetName, IMG)

S = flipud(rot90(IMG, 1));
if(islogical(S))
	S = uint8(S);
end

h5create(OutFile, DatasetName, size(S));
h5write(OutFile, DatasetName, S);
