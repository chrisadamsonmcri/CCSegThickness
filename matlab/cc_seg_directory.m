function [varargout] = cc_seg_directory(InputDir, OutputDir, varargin)

if(exist(InputDir, 'dir') ~= 7)
	error('Input directory does not exist');
end

% set up default options
DoThickness = false;
DoSeg = true;
MidSagDetect = 'acpcdetect';

z = 1;
while(z < length(varargin))
	switch(lower(varargin{z}))
		case 'midsagdetect'
			if(z == length(varargin))
				error('Argument MidSagDetect
			end
			switch
	end
end
%for z = 1:length

D = dir(fullfile(InputDir, '*.nii.gz'));

[~, ~, ~] = mkdir(OutputDir);


for z = 127:length(D)
	disp(['Subject ' num2str(z) ' of ' num2str(length(D)) ', ' D(z).name]);
	BaseName = strrep(D(z).name, '.nii.gz', '');
		cc_seg_one_subject_pipe_seg(fullfile('datasets', 'oasis_full', D(z).name), ...
 		[], ...
 		OutputDir, ...
 		[sprintf('%04d', z) '-' BaseName]);
		cc_seg_one_subject_pipe_thickness(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]));
end
delete(gcf);